module BlastParse

using FASTX: FASTA

asPercentage(s::AbstractString) = asFloat64(s) / 100.0
asFloat64(s::AbstractString) = parse(Float64, s)
asInt(s::AbstractString) = parse(Int, s, base=10)
asSubString(s::Union{String, SubString{String}}) = SubString(s)

const TYPES = Dict(
    asPercentage => Float64,
    asFloat64    => Float64,
    asInt       =>  Int,
    asSubString  => SubString{String}
)

# Not all possible fields are present
# the ones missing I just don't know how to parse
# probably not that hard, I'm just lazy
const FUNCTIONS = Dict(
    :qseqid     => asSubString,
    :qgi        => asSubString,
    :qacc       => asSubString,
    :qaccver    => asSubString,
    :qlen       => asInt,
    :sseqid     => asPercentage,
    :sallseqid  => asPercentage,
    :sgi        => asSubString,
    :sallgi     => asSubString,
    :sacc       => asSubString,
    :saccver    => asSubString,
    :sallacc    => asSubString,
    :slen       => asInt,
    :qstart     => asInt,
    :qend       => asInt,
    :sstart     => asInt,
    :send       => asInt,
    :evalue     => asFloat64,
    :bitscore   => asFloat64,
    :score      => asFloat64,
    :length     => asInt,
    :pident     => asPercentage,
    :nident     => asInt,
    :mismatch   => asInt,
    :positive   => asInt,
    :gapopen    => asInt,
    :gaps       => asInt,
    :ppos       => asPercentage,
    :staxid     => asInt,
    :ssciname   => asSubString,
    :scomname   => asSubString,
    :sblastname => asSubString,
    :sskingdom  => asSubString,
    :stitle     => asSubString,
)

"""A vector of accepted symbols that `parse` can take.
Most column names of BLAST are accepted"""
const ACCEPTED_SYMBOLS = sort!(collect(keys(FUNCTIONS)))

const DEFAULT_COLUMNS = (:qacc, :sacc, :pident, :bitscore)

"""
    gen_blastparse_code(cols::Tuple{Vararg{Symbol}}, name::Symbol)) -> Expr

Create an `Expr` that, when evaluated, defines a method `\$(name)(::IO)`.
This method parses a tab-separated BLAST output file, with the columns given
as `cols`. For valid column names, see the constant `ACCEPTED_SYMBOLS`.
Values in percentage, such as `sseqid` are scaled to [0.0:1.0].
The function returns a `Vector{<:NamedTuple{cols}}`.

# Examples
```julia
julia> @eval gen_blastparse_code((:qacc, :pident, :length), :my_parse)
my_parse (generic function with 1 method)

julia> open(my_parse, "blastout.tsv")[1:2]
Vector{NamedTuple{(:qacc, :pident, :length), Tuple{SubString{String}, Float64, Int64}}} with 2 elements:
  (qacc = "Segment1", pident = 0.81, length = 2105)
  (qacc = "Segment2", pident = 0.85, length = 2152)
```
"""
function gen_blastparse_code(cols::Tuple{Vararg{Symbol}}, name::Symbol)
    nfields = length(cols)
    functions = Function[FUNCTIONS[col] for col in cols]
    elT = Tuple{[TYPES[f] for f in functions]...}
    T = NamedTuple{cols, elT}
    body = quote end
    for i in 1:nfields
        assignment = :($(Symbol("t_", i)) = $(functions[i])(fields[$i]))
        push!(body.args, assignment)
    end
    return quote
        function $(name)(io::IO)
            fields = Vector{SubString{String}}(undef, $nfields)
            result = Vector{$T}()
            lineno = 0
            for line in eachline(io)
                lineno += 1
                #### In-place split!
                n = UInt(0)
                start = 1
                @inbounds for i in 1:ncodeunits(line)
                    if codeunit(line, i) == UInt8('\t')
                        n += 1
                        n >= $nfields && break
                        substr = SubString(line, start, i-1)
                        @inbounds fields[n] = substr
                        start = i + 1
                    end
                end
                if n+1 != $nfields
                    error("Expected $(n+1) fields, stopped at $n fields, line $lineno.")
                end
                @inbounds fields[n+1] = SubString(line, start, lastindex(line))

                ### Calculate each of the fields
                $body

                # Return namedtuple
                tup = $(Expr(:tuple, [Symbol("t_", i) for i in eachindex(functions)]...))
                namedtup = $(T)(tup)
                push!(result, namedtup)
            end
            return result
        end
    end
end

@eval $(gen_blastparse_code(DEFAULT_COLUMNS, :parse_default))

function _blastn(
    tmp_dir::AbstractString,
    query::AbstractString,
    subject::AbstractString
)
    target = tempname(tmp_dir)
    cmd = `blastn -query $query -subject $subject -outfmt "6 $(join(DEFAULT_COLUMNS, ' '))"`
    run(pipeline(cmd, stdout=target))
    rows = open(parse_default, target)
    return (rows, target)
end

function blastn(
    tmp_dir::AbstractString,
    query::Union{AbstractString, Vector{FASTA.Record}},
    subject::Union{AbstractString, Vector{FASTA.Record}},
)
    check_blastn()
    query = query isa AbstractString ? query : dump_fasta(tmp_dir, query)
    subject = subject isa AbstractString ? subject : dump_fasta(tmp_dir, subject)
    _blastn(tmp_dir, query, subject)
end

function check_blastn()
    run(pipeline(`which blastn`, stdout=devnull))
end

function dump_fasta(tmp_dir::AbstractString, recs::Vector{FASTA.Record})
    name = tempname(tmp_dir)
    open(FASTA.Writer, name) do writer
        for rec in recs
            write(writer, rec)
        end
    end
    return name
end

end # module
