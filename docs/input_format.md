# Input File Format

MitoCAPE accepts CSV files containing editing frequencies or mutation percentages by mitochondrial position.

## Required Columns

| Column | Description |
|--------|-------------|
| Position | mtDNA genomic position |
| Reference | Reference nucleotide |
| Variant | Edited nucleotide |
| Percentage | Editing percentage |

## Example

```csv
Position,Reference,Variant,Percentage
1491,C,T,32.4
1559,C,T,64.1
1561,C,T,28.7
````

## Accepted Variants

Column names are matched flexibly where possible.

## Notes

* Percentages should be numeric
* Positions should be integers
* CSV format recommended
* Header row required
