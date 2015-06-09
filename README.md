# DendroBites
This is a repository of example scripts that demonstrate how to use
[dendropy 4](http://dendropy.org/) to accomplish some basic tasks.
The goal is for this to host some small, self-contained scripts that demonstrate
best practices, but which are too arcane to belong inside DendroPy itself.


# Conventions
Code should run on python 2.7 and python 3
Each script should be (somewhat) useful in and of itself, and we should try to
limit dependencies as much as possible.

The handling of arg-parsing should be done in the `if __name__ == '__main__':`
block of code.

If there are additional user-interface steps (e.g. creation of intermediate
directories, or translation of exceptions into `sys.exit` calls) they should
be done in a `_main` function.

As much as possible, the `_main` function should delegate the work to be done
to a function with a name that corresponds with the script name.

Adhering to these conventions will make it easier for users to navigate
the examples (because of the consistency) and allow them to piece together
larger pipelines (by calling the same business logic functions that `_main`
calls). If you think that your script might be generally useful, you may want
to add it to the list of imports in `dendrobites/__init__.py` so that client
code can easily import it.

## Becoming a contributor

If you are interested in contributing scripts to `DendroBites`, feel free to email
Mark or Jeet.

# Acknowledgements

As well as depending on DendroPy, DendroBites also borrows code from it.
DendroPy is written by Jeet Sukumaran and Mark T. Holder.

