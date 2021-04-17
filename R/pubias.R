# This is the main script which combines all other functions to one package
pubias <- function(data, approach = "replication", gmm = FALSE) {

  if (gmm == FALSE){
    if (approach == "replication"){
      mle_replication()
    } else if (approach == "meta"){
      mle_meta()
    }
  } else {
    if (approach == "replication"){
      gmm_replication()
    } else if (approach == "meta"){
      gmm_meta()
    }
  }

}

# Test
# pubias <- function(data, approach = "replication", gmm = FALSE) {
#
#   if (gmm == FALSE){
#     if (approach == "replication"){
#       print("mle_replication")
#     } else if (approach == "meta"){
#       print("mle_meta")
#     }
#   } else {
#     if (approach == "replication"){
#       print("gmm_replication")
#     } else if (approach == "meta"){
#       print("gmm_meta")
#     }
#   }
#
# }
