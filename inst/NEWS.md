CHANGES IN VERSION v0.99.0
------------------------

### NEW FEATURES

* The new 'createEnrichMap()' function enables the creation of an enrichment map from enrichment results.


CHANGES IN VERSION v0.0.4
------------------------

### NEW FEATURES

* When 'intersection' column present in the GOST object, the software uses the column to create the output for Cytoscape, without having to call 'gconvert()'. The output is generated much more rapidly.


CHANGES IN VERSION v0.0.2
------------------------

### NEW FEATURES

* Added a `NEWS.md` file to track changes to the package.

### SIGNIFICANT USER-VISIBLE CHANGES

* `createNetwork()` method has a new parameter `fileName` that enables the creation of a CX JSON file when Cytoscape is not running.

### BUG FIXES

* `createCytoscapeCXJSON()` method does not replicate anymore the genes that are associated to more than one term.
