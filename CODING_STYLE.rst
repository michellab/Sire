****
`Sire <http://siremol.org>`__
****

Coding Style
============

To keep the code consistent, there a strict coding style is used in Sire.
Please follow these rules to maintain this style. Also, if you see code
that doesn't follow these rules (we are all not perfect!) then please
feel free to either correct that code, or to notify one of the lead
developers. This coding style has evolved over many years and many
developers, and is needed to ensure that the Sire source code looks
like a single homogenous set, which anyone can edit at any point.

All indentation should be multiples of 4 spaces. All tabs should be replaced by 4 spaces. 
--------------------------

For example

```c++
void foo()
{
    if (true)
    {
        for (int i=0; i<10; ++i)
        {
            //do something
        }
    }
}
```

not

```
void foo()
{
 if (true)
 {
	for (int i=0; i<10; ++i){ /* do something */ }
 }
}
```

Curly brackets should be used for all blocks, with '{' on a new line
-------------------

For example

```c++
for (int i=0; i<10; ++i)
{
    if (condition)
    {
        //do something
    }
}
```

not

```c++
for (int i=0; i<10; ++i){
    if (condition) /* do something */;
}
```


(2) Classes should be named using capital letters, using only
    the letters A-Za-z and numbers 0-9. Please do not use underscores.
    For example BigMolecule is acceptable, but bigMolecule, Bigmolecule
    Big_Molecule or bigmolecule are not.

(3) Functions (methods) should be named in the same way, except that
    the first letter should not be capitalised. For example
    getRadius() is acceptable, but GetRadius(), getradius() or
    get_radius() is not.

(4) Variables (member data) should be named using all small case
    letters or numbers. Underscores should be used to separate
    words, and obvious abbreviations are recommended (e.g. mol for molecule).
    For example, added_mol is acceptable, but added_molecule should be avoided,
    and Added_Mol, addedMol, Addedmol are all not acceptable

(5) Exceptions are named in the same way as variables, except abbreviations
    should not be used, e.g. missing_molecule is acceptable, but missing_mol is not

(6) No line should be over 90 characters long. Long lines should be split,
    with the extra part indented so that it lines up with the above line, e.g.

    AtomCoords coords = mol.atom( AtomName("O00") )
                           .property("coordinates")
                           .asA<AtomCoords>();

(7) Which brings me on to - always code using a fixed-width font. The code
    uses whitespace and indentation to make things clear, and this is lost
    if you use a variable width font

(8) Use whitespace to make the code clean. For example, always have a blank
    line before a code block (e.g. function, if statement, for loop),
    except if it comes directly after an open brace "{". For example

void foo()
{
    int a;

    if (a == 5)
    {
        for (int i=0; i<10; ++i)
        {
            if (b == 10)
            {
                a = 5 * b;

                for (int j=0; j<11; ++j)
                {}
            }
        }
    }
}

(9) Speaking of braces, please use the above style - e.g. braces are on their own
    line and line up. This makes it much easier to read.

(10) Sire uses doxygen to autogenerate the API documentation. This means comments
     should be written using these rules;
  
     (i) All class and function comments should start /** and end with */
     (ii) If you author a class, add a @author Your Name to the class comment
     (iii) If you function throws an exception, add a \throw Namespace::exception 
           to the function comment
     (iv) Use "//" for all other comments (even multiline). This is so that it is
          possible to quickly comment out blocks of text using "/*" and "*/"
