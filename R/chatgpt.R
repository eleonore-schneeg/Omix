#' Calls chatGPT with a specific prompt
#'
#' @param prompt Character vector
#' @param api_key Personal API key
#'
#' @return ChatGPT answer
#' @export
#' @importFrom httr add_headers
#' @importFrom stringr str_trim

ask_chatgpt <- function(prompt,
                        api_key="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") {
  response <- POST(
    url = "https://api.openai.com/v1/chat/completions",
    httr::add_headers(Authorization = paste("Bearer", api_key)),
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
  stringr::str_trim(content(response)$choices[[1]]$message$content)
}
