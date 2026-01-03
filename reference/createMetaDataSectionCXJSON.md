# Create meta data section for the CX JSON file

Create meta data section for the CX JSON file that contains the network
information

## Usage

``` r
createMetaDataSectionCXJSON()
```

## Value

a `JSON` object that contains the meta data section related to the
network

## Author

Astrid Deschênes

## Examples

``` r

## Create the JSON object that contains the meta data information
enrichViewNet:::createMetaDataSectionCXJSON()
#> {"metaData":[{"name":"nodes","version":"1.0"},{"name":"edges","version":"1.0"},{"name":"edgeAttributes","version":"1.0"},{"name":"nodeAttributes","version":"1.0"},{"name":"cyHiddenAttributes","version":"1.0"},{"name":"cyNetworkRelations","version":"1.0"},{"name":"cyGroups","version":"1.0"},{"name":"networkAttributes","version":"1.0"},{"name":"cyTableColumn","version":"1.0"},{"name":"cySubNetworks","version":"1.0"}]} 
```
