# math_group_tools
Basic tools for working with groups, (abstract algebra)

## Example Input and Features:
For help and info:

`python groups.py -h`

`python groups.py -i`

`python groups.py -i -g GL`

### Easy defining of groups:

`python groups.py -g MatrixGroup "[[[0,0]],[[0,1]],[[0,2]],[[0,3]],[[1,0]],[[1,1]],[[1,2]],[[1,3]],[[2,0]],[[2,1]],[[2,2]],[[2,3]],[[3,0]],[[3,1]],[[3,2]],[[3,3]]],op=<matrixelement,<addmod,4>>,name=Z4xZ4"`

### Control types and ordering of tasks:

`python groups.py -g S 3 -t abelian orders "centralizer{S,3}" "lcosets{PermutationGroup,[[1,2],[]]}" cache "rcosets{PermutationGeneratorGroup,[[1,2]]}" cache`

### Multiple groups tested at once:

`python groups.py -g S 4 -g M "2,{Z,2},<matrixelement,<addmod,2>,cache=256>"`

`python groups.py -g PermutationGeneratorGroup "[[[1,2,3]],[[1,2]]]" -g S 3`
