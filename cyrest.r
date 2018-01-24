### functions to start up cyrest and its functions and establish a connnection with Cytoscape
## source this file but make sure Cytoscape is ready to use before that

## load libraries
library(RJSONIO)
library(igraph)
library(httr)

## functions necessary for cyrest
source('toCytoscape.R')

## default port for Cytoscape
port.number = 1234
base.url = paste("http://localhost:", toString(port.number), "/v1", sep="")
version.url = paste(base.url, "version", sep="/")
cytoscape.version = GET(version.url)
cy.version = fromJSON(rawToChar(cytoscape.version$content))
network.url = paste(base.url, "networks", sep="/")
## Seeing this indicates that the connection is ready to use, so only need to do once
print(cy.version)