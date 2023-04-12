#' Calls chatGPT with a specific prompt
#'
#' @param prompt Character vector
#' @param api_key Personal API key
#'
#' @return ChatGPT answer
#' @family Helper
#' @export
#' @importFrom httr add_headers
#' @importFrom stringr str_trim
ask_chatgpt <- function(prompt='What is multi-omics',
                        api_key=NULL) {
  if(api_key==NULL){
    "No OpenAI API key provided, please fill api_key parameter to proceed"
  }
  
  if(!is.null(api_key)){
  response <- POST(
    url = "https://api.openai.com/v1/chat/completions", 
    add_headers(Authorization = paste("Bearer", api_key)),
    content_type_json(),
    encode = "json",
    body = list(
      model = "gpt-3.5-turbo",
      messages = list(list(
        role = "user", 
        content = prompt
      ))
    )
  )
  ret <- content(response)
  if ("error" %in% names(ret)) warning(ret$error$message)
  if ("message" %in% names(ret$choices[[1]]))
    cat(stringr::str_trim(ret$choices[[1]]$message$content))
  return(invisible(ret))
  }
  
}

