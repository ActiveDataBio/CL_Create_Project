# get arguments (rdata file path)
args <- commandArgs(trailingOnly = TRUE)
#print(args)
setwd(args[1])
load(args[2])
hm <- hc.out

hc_row <- as.hclust( hm$rowDendrogram )
hc_col <- as.hclust( hm$colDendrogram )

library(rjson)

#convert output from hclust into a nested JSON file
HCtoJSON<-function(hc) {
  labels<-hc$labels
  merge<-data.frame(hc$merge)
  height<-as.vector(hc$height)
  
  for (i in (1:nrow(merge))) {
    if (merge[i,1]<0 & merge[i,2]<0) {
      eval(parse(text=paste0("node", i, "<-list(name=\"node", i, "\", children=list(list(name=labels[-merge[i,1]]),list(name=labels[-merge[i,2]])))")))
    } else if (merge[i,1]>0 & merge[i,2]<0) {
      if (i > 1 && height[i] - height[merge[i,1]] < 0.0001) {
        eval(parse(text=paste0("node", i, "<-list(name=\"node", i, "\", children=list(list(name=labels[-merge[i,2]])))")))
        eval(parse(text=paste0("numChild <- length(node", merge[i,1],"$children)")))
        for (j in (1:numChild)) {
          eval(parse(text=paste0("node", i, "$children[j+1]<-node", merge[i,1], "$children[j]")))
        }
      } else {
        eval(parse(text=paste0("node", i, "<-list(name=\"node", i, "\", children=list(node", merge[i,1], ", list(name=labels[-merge[i,2]])))")))
      }
    } else if (merge[i,1]<0 & merge[i,2]>0) {
      if (i > 1 && height[i] - height[merge[i,2]] < 0.0001) {
        eval(parse(text=paste0("node", i, "<-list(name=\"node", i, "\", children=list(list(name=labels[-merge[i,1]])))")))
        eval(parse(text=paste0("numChild <- length(node", merge[i,2],"$children)")))
        for (j in (1:numChild)) {
          eval(parse(text=paste0("node", i, "$children[j+1]<-node", merge[i,2], "$children[j]")))
        }
      } else {
        eval(parse(text=paste0("node", i, "<-list(name=\"node", i, "\", children=list(list(name=labels[-merge[i,1]]), node", merge[i,2],"))")))
      }
    } else if (merge[i,1]>0 & merge[i,2]>0) {
      if (i > 1 && height[i] - height[merge[i,1]] < 0.0001) {
        eval(parse(text=paste0("node", i, "<-list(name=\"node", i, "\", children=list())")))
        eval(parse(text=paste0("node", i, "$children<-node", merge[i,1], "$children")))
        eval(parse(text=paste0("num1 <- length(node", i,"$children)")))
        eval(parse(text=paste0("num2 <- length(node", merge[i,2],"$children)")))
        for (j in (1:num2)) {
          eval(parse(text=paste0("node", i, "$children[j+num1]<-node", merge[i,2], "$children[j]")))
        }
      } else {
        eval(parse(text=paste0("node", i, "<-list(name=\"node", i, "\", children=list(node",merge[i,1] , ", node" , merge[i,2]," ))")))
      }
    }
  }
  
  eval(parse(text=paste0("JSON<-toJSON(node",nrow(merge), ")")))
  return(JSON)
}

JSON_row<-HCtoJSON(hc_row)
JSON_col<-HCtoJSON(hc_col)

fileConn<-file("dendro_row.json")
writeLines(paste0(JSON_row), fileConn)
close(fileConn)

fileConn<-file("dendro_col.json")
writeLines(paste0(JSON_col), fileConn)
close(fileConn)

write.csv(file="matrix.csv", x=hm$carpet)