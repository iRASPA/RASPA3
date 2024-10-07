Largely based on Google C++ Style Guide:
(https://sun.iwu.edu/~mliffito/cs_codex/posts/google-c++-style-guide/)

Naming
======
1) Local and global variables will be named in lower_case_with_underscores.
2) Structure member variables (attributes) will be named with_a_trailing_underscore_.
3) Structure names will be named in upper CamelCase
4) Functions will be named in upper CamelCase
5) Enums, will be named in upper CamelCase
6) Structure methods will be named in lower CamelCase
7) All names should be descriptive, other than loop counter variables like i. 

Comments
========
1) Variable, class, and function comments come before the commented component.
2) Doxygen comments:
/// Inserts a molecule into the vector of atoms.
///
/// Note: updates the numberOfMoleculesPerComponent, numberOfIntegerMoleculesPerComponent,
///       numberOfPseudoAtoms, totalNumberOfPseudoAtoms.
/// - Parameters:
///   - selectedComponent: the index of the component
///   - atoms: vector of atoms to be inserted
/// - returns:
3) Prefer // over /* and */. 

Formatting
==========
1) Terminal: 120 columns
2) Braces: Allman style
3) Horizontal whitespace: Use 1 space around (most) operators. E.g., x = y + z as opposed to x=y+z. Exceptions: ., ->, ++, and --. 
4) Do not put a space between a function name and its argument or parameter list. E.g., write Function(arg1, arg2), not Function (arg1, arg2).
5) Vertical whitespace: Use 2 blank lines between functions/classes etc., as well as between related sections of code within a function (sort of like paragraph breaks in writing).
6) Indentation: Indent each nested block 2 spaces at a time. Do not use tabs in your code. You should set your editor to emit spaces when you hit the tab key.

Writing structures
==================
1) Structures should be declared and implemented in their own files, which are named structure_name.ixx (declaration) and structure_name.cpp (implementation).
2) Structures data members (attributes) are public
3) All structures should have explicit constructors and destructors, even if they do nothing, for clarity.
4) Use this-> to access an object’s attributes and methods from within its methods — this immediately differentiates them from normal variables and functions.

Type deduction
==============
1) avoid 'auto' to prioritize code clarity and safety over convenience.

Definitions
===========
upper camel case, or Pascal case: ThisIsUpperCamelCase or ThisIsPascalCase
lower camel case: thisIsLowerCamelCase
snake case: this_is_snake_case

Allman style: This style puts the brace associated with a control statement on the next line, indented to the same level as the control statement. Statements within the braces are indented to the next level.
while (x == y)
{
  foo();
  bar();
}


