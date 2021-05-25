
# gprofiler2cytoscape v0.0.2

NEW FEATURES

* Added a `NEWS.md` file to track changes to the package.

SIGNIFICANT USER-VISIBLE CHANGES

* `createNetwork()` method has a new parameter `fileName` that enables the creation of a CX JSON file when Cytoscape is not running.

BUG FIXES

* `createCytoscapeCXJSON()` method does not replicate anymore the genes that are associated to more than one term.