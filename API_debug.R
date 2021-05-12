library(httr)
library(jsonlite)

genes <- read_csv("cluster7_previous_RNA.csv") %>%
  pull(symbol) %>%
  distinct()

res <- POST(url = "https://reactome.org/AnalysisService/identifiers/projection/?interactors=false&species=Homo%20sapiens&pageSize=200&page=1&sortBy=ENTITIES_FDR&order=ASC&resource=TOTAL&pValue=1&includeDisease=true",
            body = genes,
            accept_json(),
            content_type("text/plain"))

raw_content <- res$content
char_from_raw <- rawToChar(res$content)
content_json_char_from_raw <- fromJSON(rawToChar(res$content))
pathways <- fromJSON(rawToChar(res$content))$pathways
