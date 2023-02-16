
1. 3s_6loci_fixnet.xml

This BEASTInterface (initEmbed) has no input with name taxonset. Choose one of these inputs: speciesNetwork,geneTree,operator,weight

Error detected about here:

```xml
<init id='SNI' spec='speciesnetwork.SpeciesNetworkInitializer'>
   <rebuildEmbedding id='initEmbed' spec='speciesnetwork.operators.RebuildEmbedding'>
```

Attempting to fix this XML, but it looks `speciesnetwork.GeneTreeInSpeciesNetwork` does not exist: https://github.com/zhangchicool/speciesnetwork/tree/179b241d47edba679c3b975a17a2f9279a1039d0/src/speciesnetwork
```xml
<geneTreeWithin id="geneTree.t:locus4" spec="speciesnetwork.GeneTreeInSpeciesNetwork" geneTree="@Tree.t:locus4" speciesNetwork="@Network.t:Species"/>
```


2. 3s_6loci_infernet.xml

This BEASTInterface (initEmbed) has no input with name taxonset. Choose one of these inputs: speciesNetwork,geneTree,operator,weight

Error detected about here:
```xml
<init id='SNI' spec='speciesnetwork.SpeciesNetworkInitializer'>
   <rebuildEmbedding id='initEmbed' spec='speciesnetwork.operators.RebuildEmbedding'>
```

Due to https://github.com/zhangchicool/speciesnetwork/commit/c25a796a915ddc62bd527da1c1beb8178c09cf53#diff-da3bbecb13258d35dcc77d2c16d1a76d8cbeffd5aa7cc90b6da0cb9d80680f6d

