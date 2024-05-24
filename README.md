# BlastParse.jl
Lazy parsing of tabular BLAST files.

## Usage:
* Parse default tabular BLAST output:
```
using BlastParse: each_hit

open("blast_output.txt") do io
    for hit in each_hit(io)
        ...
    end
end
```

* Note that fractions that are written as percent in the BLAST output, like identity,
  are always parsed to a float in the correct 0-1 range.
* To parse BLAST with non-default columns, see the docstring of `each_hit`
