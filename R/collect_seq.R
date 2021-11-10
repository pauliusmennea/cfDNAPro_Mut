
#' collect_sequence
#'
#' @param method 
#' @param frags 
#' @param pos 
#' @param offset 
#'
#' @return
#' @export
#'
#' @examples
collect_sequence <- 
  function(method,
           frags,
           pos,
           offset){
    if(method=="flank"){
      if(pos == "start"){
        collectedseq=getSeq(genome,flank(unstrand(frags),offset+1,start=TRUE))
      }
       if(pos == "end"){
        collectedseq=getSeq(genome,flank(unstrand(frags),offset+1,start=FALSE))
      }     
    }
    if(method=="resize"){
      if(pos == "start"){
        collectedseq=getSeq(genome,resize(unstrand(frags),offset,fix="start"))
      }
      if(pos == "end"){
        collectedseq=getSeq(genome,resize(unstrand(frags),offset,fix="end"))
      }      
    }
    return(collectedseq)
  }

