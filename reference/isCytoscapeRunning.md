# Verifying that Cytoscape is running

Verifying that Cytoscape is running

## Usage

``` r
isCytoscapeRunning()
```

## Value

a `logical` indicating if Cytoscape is running.

## Author

Astrid Deschênes

## Examples

``` r

## Test if Cytoscape is running
enrichViewNet:::isCytoscapeRunning()
#> Unable to connect to Cytoscape. 
#> CX JSON file will be created.
#> [1] FALSE
```
