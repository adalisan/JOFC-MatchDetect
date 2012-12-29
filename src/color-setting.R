# Example preprocessing script.

 require(RColorBrewer)

	colors.vec <- c("red","green","aquamarine","purple",
			"darkblue","salmon","rosybrown","magenta","orange")
	colors.vec[3]<-"gold4"
	colors.vec[2]<-"darkblue"
	colors.vec[4]<-"darkorange4"
	colors.vec[9]<-"red"
	colors.vec.len<-length(colors.vec)
	
      colors.vec.len <- w.val.len

	colors.vec<- brewer.pal(colors.vec.len,"Spectral")
    #colors.vec.len <- w.val.len+3
     # colors.vec <- colors.vec[-(1:3)]
      colors.vec.len <- length(colors.vec)
      
	
	colors.vec[colors.vec.len+1]<-"cornflowerblue"
	colors.vec[colors.vec.len+2]<-"red"
	colors.vec[colors.vec.len+3]<-"cyan"
	
	
	
	
	
	colors.vec.len <- length(colors.vec)
	colors.vec[(colors.vec.len-2):colors.vec.len] <- brewer.pal(3,"Set3")
	#palette("YlOrRd")