# Copyright (C) 2018 Christopher Cave-Ayland
# You may use, distribute and modify this code under the terms of the
# MIT open source license. See https://opensource.org/licenses/MIT if
# a copy was not provided with this code.

array_help_string = """This object provides a pythonic interface to an underlying C++ array,
mimicing the semantics of a list. This object has methods and attributes that
wrap calls to the openmm base class, {0}.
This object can be indexed to call the get and set routines detailed below.
Where possible, accessed array elements will be provided as a namedtuple
to give additional information about the underlying data.
This function may also have append and pop attributes to add or remove
elements of the array.

This object provides a list-like interface that wraps the below
function calls from the openmm base class:
Wrapped length function - {1}
Wrapped get functions - {2}
Wrapped set functions - {3}
Wrapped add function  - {4}
Wrapped remover function - {5}
"""

append_help_string = """This function acts as a wrapper to a method that adds an element to a
C++ array. This function expects a single tuple of values that are passed as
arguments to the function {0}.

The help string of {0} is provided below, detailing
the values that should be given in the tuple passed to this function -

{1}
"""

pop_help_string = """This function acts as a wrapper to a method that removes an element
from an underlying C++ array. This function expects an integer index value
that is passed to the function {0}.
Unlike list.pop this value does not return the removed element.

The help string of {0} is provided below, detailing
the effect of removing an element -

{1}
"""

dict_help_string = """hehehehehe"""

class_help_string = """This class provides a pythonic wrapper for the openmm base class
 {0}. 

Wherever possible 'get', 'set' and 'add' methods from the base class have
been combined into a single interface that provides list-like semantics.

The help string of {0} is provided below -

{1}
"""
