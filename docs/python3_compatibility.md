# Python 3 compatibility

## Using the PyChar replace with regex

find (with regex ticked) : ([a-z,A-Z,.]*).has_key\(\'([a-z,A-Z,_,0-9]*)\'\)
replace by : '$2' in $1

find (with regex ticked) : not '([a-z,A-Z,_,0-9]*)' in self.data
replace by : '$1' not in self.data