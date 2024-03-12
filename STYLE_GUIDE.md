Terminal: 120 columns
Braces: Allman style
Tabs Versus spaces: Use only spaces, and indent 2 spaces at a time
File Names: Filenames should be all lowercase and can include underscores
Extensions: source.cpp, source.ixx
Efficiency: Use ++i as opposed to i++. (https://medium.com/better-programming/stop-using-i-in-your-loops-1f906520d548)

Integers: <cstdint> defines types like int16_t, uint32_t, int64_t, etc. You should always use those in 
          preference to short, unsigned long long and the like, when you need a guarantee on the size of an integer.

Type deducation: avoid 'auto' to prioritize code clarity and safety over convenience.
Exception: avoid, use std::expected
Annotation: [[likely]], [[unlikely]], [[no_discard]], [[maybe_unused]]

Local variable names: 
The names of variables (including function parameters) and data members are all lowercase, with underscores between words. 

Structures and enums names: UpperCamelCase, should start with a capital letter and have a capital letter for each new word.
structure member names: lowerCamelCase, should start with a lowercase letter and have a capital letter for each new word.
Function calls: lowerCamelCase, should start with a lowercase letter and have a capital letter for each new word.

Function calls: If the arguments do not all fit on one line, they should be broken up into multiple lines, 
                with each subsequent line aligned with the first argument. Do not add spaces after the open 
                paren or before the close paren.
