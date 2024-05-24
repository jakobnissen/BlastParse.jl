"""
    BlastParse

Lazy parsing of tabular BLAST output files.

This module implements the convenience function `each_hit`, which is a lazy
iterator of rows of a tabular BLAST output file.

Example usage:
```julia
open("blast_output.tsv") do io
    for hit in each_hit(io)
        # do something here
    end
end
```
"""
module BlastParse

export each_hit

asPercentage(s::AbstractString) = asFloat64(s) / 100.0
asFloat64(s::AbstractString) = parse(Float64, s)
asInt(s::AbstractString) = parse(Int, s; base=10)
asIntPerc(s::AbstractString) = parse(Int, s; base=10) / 100.0
asSubString(s::Union{String, SubString{String}}) = SubString(s)

const TYPES = Dict(
    asPercentage => Float64,
    asFloat64 => Float64,
    asInt => Int,
    asSubString => SubString{String},
    asIntPerc => Float64,
)

# Not all possible fields are present
# the ones missing I just don't know how to parse
# probably not that hard, I'm just lazy
const FUNCTIONS = Dict(
    :qseqid => asSubString,
    :qgi => asSubString,
    :qacc => asSubString,
    :qaccver => asSubString,
    :qlen => asInt,
    :sseqid => asPercentage,
    :sallseqid => asPercentage,
    :sgi => asSubString,
    :sallgi => asSubString,
    :sacc => asSubString,
    :saccver => asSubString,
    :sallacc => asSubString,
    :slen => asInt,
    :qstart => asInt,
    :qend => asInt,
    :sstart => asInt,
    :send => asInt,
    :evalue => asFloat64,
    :bitscore => asFloat64,
    :score => asFloat64,
    :length => asInt,
    :pident => asPercentage,
    :nident => asInt,
    :mismatch => asInt,
    :positive => asInt,
    :gapopen => asInt,
    :gaps => asInt,
    :ppos => asPercentage,
    :staxid => asInt,
    :ssciname => asSubString,
    :scomname => asSubString,
    :sblastname => asSubString,
    :sskingdom => asSubString,
    :stitle => asSubString,
    :qcovs => asIntPerc,
    :qcovhsp => asIntPerc,
    :qcovus => asIntPerc,
)

"""A vector of accepted symbols that `parse` can take.
Most column names of BLAST are accepted"""
const ACCEPTED_SYMBOLS = sort!(collect(keys(FUNCTIONS)))

# Default fields in NCBI blastn 2.6.0
const DEFAULT_COLUMNS = (
    :qaccver,
    :saccver,
    :pident,
    :length,
    :mismatch,
    :gapopen,
    :qstart,
    :qend,
    :sstart,
    :send,
    :evalue,
    :bitscore,
)

struct HitIterator{T, I}
    lines::I
    substrings::Vector{SubString{String}}
end

each_hit(io::IO) = each_hit(io, Val{DEFAULT_COLUMNS}())

Base.@assume_effects :foldable function check_columns(cols)
    for col in cols
        if !haskey(FUNCTIONS, col)
            error("Not column name \"$(repr(col))\", but no such BLAST column type exist")
        end
    end
    nothing
end

function each_hit(io::IO, ::Val{C}) where C
    @inline check_columns(C)
    substrings = fill(SubString("", 1, 0), length(C))
    lines = eachline(io)
    HitIterator{C, typeof(lines)}(lines, substrings)
end

"""
    each_hit(io::IO, [columns::Val{<:Tuple{Vararg{Symbol}}}])

Produce a lazy iterator over each row of the tabular BLAST output
stored in `io`.

The optional `columns` argument controls the expected columns of `io`,
if not passed, the default BLAST columns will be assumed.
The argument is passed as a `Val` of tuple of symbols, e.g.
`Val{(:qaccver, :saccver, :pident, :bitscore)}()`.
For all accepted column names, see `BlastParse.ACCEPTED_SYMBOLS`.

All fractions are parsed to floats in [0, 1], e.g. a `pident` of
95.1 is parsed to the float 0.951.
"""
each_hit

@generated function Base.eltype(::Type{<:HitIterator{T}}) where T
    functions = Function[FUNCTIONS[col] for col in T]
    elT = Tuple{[TYPES[f] for f in functions]...}
    :(NamedTuple{$(T), $(elT)})
end

Base.IteratorSize(::Type{<:HitIterator}) = Base.SizeUnknown()

@generated function Base.iterate(it::HitIterator{T}, lineno::Int=1) where T
    nfields = length(T)
    functions = Function[FUNCTIONS[col] for col in T]
    body = quote end
    for i in 1:nfields
        assignment = :($(Symbol("t_", i)) = $(functions[i])(@inbounds fields[$i]))
        push!(body.args, assignment)
    end
    quote
        nfields = $nfields
        maybe_line = iterate(it.lines)
        isnothing(maybe_line) && return nothing
        (line, _) = maybe_line
        fields = it.substrings
        n = 0
        start = 1
        @inbounds for i in 1:ncodeunits(line)
            if codeunit(line, i) == UInt8('\t')
                n += 1
                n >= $nfields && break
                substr = SubString(line, start, prevind(line, i))
                fields[n] = substr
                start = i + 1
            end
        end
        if n + 1 != $nfields
            error(lazy"On line $lineno, expected $(nfields) fields, stopped at $(n + 1) fields")
        end
        @inbounds fields[n + 1] = SubString(line, start, lastindex(line))

        ### Calculate each of the fields
        $body

        # Return namedtuple
        tup = $(Expr(:tuple, [Symbol("t_", i) for i in eachindex(functions)]...))
        namedtup = eltype(typeof(it))(tup)
        (namedtup, lineno + 1)
    end
end

end # module
