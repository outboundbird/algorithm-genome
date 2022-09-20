FindMyNeighbor.2 <- function(data){
  
  chr= data.frame(data, stringsAsFactors = F)
  
  for(col in 1:2){
    ifelse(is.numeric(chr[,col]) , names(chr)[col] <-"pos", names(chr)[col] <-"source")
  }  
  
  
  # sort data first on position then on source - maker sure 'inf' before 'ref' 
  chr = chr[order(chr$pos, chr$source),]
  rownames(chr)<-NULL
  idx_mark= as.numeric(rownames(chr[chr$source=="inf",]))
  i = idx_mark[1] 
  
  # define when the loop will stop: when the index hit the last position on the list. Inside the loop the value of i keeps updating until it run out the last position. the reason didn't set as the last position is that I want to include the last position in the output  
  last.pos <- nrow(chr) 
  this_is_the_end <- i > last.pos
  
  rest <- list()
  # ==================== begin of the while loop =================
  
  while(!this_is_the_end) {
    
    # check if this marker is inference or reference. it should be cuz i define the i to be the 1st pos of inference. And it should be consistent with the end!!!<======== (not necessarily if at i+1 row is a reference) 
    is.marker <- chr[i,"source"] == "inf"
    # if marker is reference, find the next inference marker 
    if(!is.marker){
      
      warning(paste("Row",i,"is not a inference marker."))
      
      while(!is.marker & i<= last.pos){
        is.marker <- chr[i, "source"] =="inf"
        if(is.marker ) break
        i = i + 1
      }
    }
    
    if(i < last.pos){
      # i should indiate the row of inference marker now 
      # define the starting point of the inference section i-1 is the backward neighbor reference
      start_m = i
      
      k = i
      # check if the next position of i is inference marker, if FALSE, stop; if TRUE, continue to the next row
      
      while(is.marker & k < last.pos){
        k = k +1
        is.marker <- chr[k,"source"]=="inf"
        if(k == last.pos) break
      }
      
      # this define the forward end of the inference section 
      if(k == last.pos & is.marker){
        stop_m = k
      }else{
        stop_m = k-1
      }
      
      # find the backward and forward neighbors
      # backward reference neighbor
      if(i==1){
        ref_b_pos= 0
      }else if(i > 1 & chr[(i-1), "source"]=="ref") {
        ref_b_pos = chr[(i-1), "pos"]
      }else{
        stop(paste("Can't find the backward reference neighbor. Check row number:", i))
      }
      
      # forward reference neighbor
      if(is.marker & k == last.pos){
        warning(paste("No last forward reference neighbor; end of the list:", k,". Set forward neighbor to 0."))
        ref_f_pos = 0
      } else if(!is.marker){
        ref_f_pos = chr[k,"pos"]
      }else{
        stop(paste("Error check row number:", k,". Not reference neighbor."))
      }
      
      # to avoid zero values (inf marker falls on the ref markers)
      m = k + 1 ; # test if the next position is a ref to make k +1 
      
      if(m < last.pos){
        is.marker <- chr[m, "source"]=="inf"
        while(is.marker){
          m = m +1
          is.marker <- chr[m,"source"]=="inf"
        }
        ref_ff_pos = chr[m, "pos"]; # forward to the next 'ref' row to avoid 0 value
      }else{
        ref_ff_pos = NA
      }
      
      # calculate the distance between the two neighbors and find the closest neighbor 
      mark_pos = chr[start_m:stop_m,"pos"] 
      ll = abs(ref_b_pos - mark_pos)
      ul = abs(ref_f_pos - mark_pos)
      ufl = abs(ref_ff_pos - mark_pos); # upper 'forward' limit
      dist= cbind(ll, ul, ufl)
      
      output <- sapply(1: nrow(dist), function(row){
        # find non-zero minimum distance to the neighbor reference
        # if the inf marker falls on the ref marker, this step will find the next nearest ref marker
        r.min <- min(dist[row,dist[row,] >0], na.rm =T)
        idx.min <-match(r.min, dist[row,])
        site.close <-c(ref_b_pos, ref_f_pos, ref_ff_pos)[idx.min]
        
        site.dist <- min(dist[row,][dist[row,]>0], na.rm=T)
        mark <- mark_pos[row]
        return(c(mark,site.dist, site.close))
      })
      
      rest[[length(rest)+1]] <- t(output)
      
      # stop signal i = last position in the list
      # if i is at the last row of the data
      i= k+1
      this_is_the_end <- i > last.pos
    }else{
      this_is_the_end <- i > last.pos
    }
  }
  
  # ================== end of the loop =======================
  map = do.call('rbind', rest)
  return(map)
}
