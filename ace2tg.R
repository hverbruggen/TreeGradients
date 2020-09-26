ace2tg <- function(ace, phy, data, file = "") {
	
	library(geiger)
	
	if (class(ace) != "ace") print("error -- first argument is not of class ace")
	if (class(phy) != "phylo") print("error -- second argument is not of class phylo")
	
	cat("",file = file)
	
	for (name in names(data)) {
		cat(paste(data[name],"\t",name,"\n", sep = ""), file = file, append = TRUE)
	}

	for (node in names(ace$ace)) {
		cat(paste(ace$ace[node],"\t",sep = ""), file = file, append = TRUE)
		cat(paste(node.leaves(phy,as.integer(node)),sep = " "), file = file, append = TRUE)
		cat("\n", file = file, append = TRUE)
	}
	
}
