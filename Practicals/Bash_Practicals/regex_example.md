# Some correct options for the regular expression in getLocations.sh  #

1. `sed -r 's/.+BDGP6:([^:]*):([0-9]+):([0-9]+).+gene:([^ ]+).+/\4\t\1\t\2\t\3/g' `

or 

2. `sed -r 's/.+BDGP6:([^:]+):([0-9]+):([0-9]+).+gene:([^ ]+).+/\4\t\1\t\2\t\3/g' `

or 

3. `sed -r 's/.+BDGP6:(.*):([0-9]+):([0-9]+):.+gene:([^ ]+).+/\4\t\1\t\2\t\3/g' `

or

4.  `sed -r 's/.+BDGP6:(.+):([0-9]+):([0-9]+):.+gene:([^ ]+).+/\4\t\1\t\2\t\3/g' `


## An  incorrect option ##

5.  `sed -r 's/.+BDGP6:(.*):([0-9]+):([0-9]+).+gene:([^ ]+).+/\4\t\1\t\2\t\3/g' `

This fails because the first `.*` is greedy and keeps going past the first `:`. This gets things out of step with the data. This doesn't happen with `:` before `.+`. 

To see this you have to know that `:(.+):([0-9]+):([0-9]+)` matches ":2L:21215306:21215178", but it also matches ":2L:21215306:21215178:1", by making the first capture "2L:21215306". **This is only possible in the case that the last colon-separated field is digits only**. When it has a "-", this is not true, so only the first option is available.

Alternatively, you can make it impossible for the first match to overrun by using `([^:]*)` as in the first two examples. 