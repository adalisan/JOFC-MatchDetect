# Example preprocessing script.



	colors.vec <- c("red","green","aquamarine","purple",
			"darkblue","salmon","rosybrown","magenta","orange")
	colors.vec[3]<-"gold4"
	colors.vec[2]<-"darkblue"
	colors.vec[4]<-"darkorange4"
	colors.vec[9]<-"red"
	colors.vec.len<-length(colors.vec)
	
  require(RColorBrewer)
	colors.vec<- brewer.pal(colors.vec.len,"YlOrRd")
	
	colors.vec[colors.vec.len+1]<-"cornflowerblue"
	colors.vec[colors.vec.len+2]<-"azure3"
	colors.vec[colors.vec.len+3]<-"cyan"
	
	
	
	
	
	colors.vec.len <- length(colors.vec)
	colors.vec[(colors.vec.len-2):colors.vec.len] <- brewer.pal(3,"Set3")
	#palette("YlOrRd")