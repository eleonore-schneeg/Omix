OpenTargets = function(disease_id,size= 2000){
  query_url = 'https://api.platform.opentargets.org/api/v4/graphql'
  
  # Building query:
  request_body = list(
    operationName= 'DiseaseAssociationsQuery',
    variables = list(
      efoId= disease_id,
      index= 0,
      size= size,
      sortBy= 'score',
      filter= '',
      aggregationFilters = list()
      
    ),
    query = '
            query DiseaseAssociationsQuery($efoId: String!, $index: Int!, $size: Int!, $filter: String, $sortBy: String!, $aggregationFilters: [AggregationFilter!]) {
                disease(efoId: $efoId) {
                    name
                    associatedTargets(page: {index: $index, size: $size}, orderByScore: $sortBy, BFilter: $filter, aggregationFilters: $aggregationFilters) {
                        count
                        rows {
                            target{
                                id
                                approvedSymbol
                            }
                            score
                            datatypeScores{
                                id
                                score
                            }
                        }
                    }
                }
            }
        '
  )
  # Retrieve data:
  response = httr::POST(url=query_url, body=request_body, encode='json')
  
  # Parse data:
  char = rawToChar(response$content)
  data = jsonlite::fromJSON(char)
  
  # Extracting associations:
  associations = data$data$disease$associatedTargets$rows
  
  # Adding target and disease columns to dataframe:
  associations$diseaseName = data$data$disease$name
  associations$targetId  = associations$target$id
  associations$targetSymbol = associations$target$approvedSymbol
  
  # Dropping unused columns and return:
  return (associations[, c('targetId', 'targetSymbol', 'diseaseName', 'score', 'datatypeScores')])
}


get_diseaseAssociations_df<-function(disease_id = 'MONDO_0004975', size=1000){
  
  library(dplyr)
  library(ghql)
  library(jsonlite)
  library(reshape2)
  
  data=OpenTargets(disease_id,size=size)
  data=as.data.frame(data)
  data=as.data.frame(data)%>% 
    flatten()%>%
    rename(overallScore = score) %>% 
    tidyr::unnest(datatypeScores) %>% 
    tidyr::pivot_wider(names_from = "id", values_from = "score")
  data=reshape2::melt(data)
  return(data)
}


#data=get_diseaseAssociations_df(disease_id = 'MONDO_0004975')

filter_OpenTarget<-function(opentarget_results=data,genes){
  library(reshape2)
  final=opentarget_results[which(opentarget_results$targetSymbol %in% genes),]
  return(final)
}

plot_filter_OpenTarget<-function(opentarget_results=data,genes){
  library(reshape2)
  library(plotly)

  
  final=opentarget_results[which(opentarget_results$targetSymbol %in% genes),]
  # Heatmap 
  p=ggplot(final, aes(variable, targetSymbol, fill= value)) + 
    geom_tile()+
    labs(x="Association scores from OpenTargets", y=paste("Targets associated with",(final$diseaseName)), title="") +
    scale_fill_gradient(low="white",high="darkblue")+ 
    guides(fill=guide_colorbar("Score (%)")) +
    theme(
      axis.text.x = element_text(size=rel(1.0), color="black"),
      axis.text.y = element_text(size=rel(1.1), color="black"),
      axis.line = element_blank(),
      axis.ticks =  element_blank(),
      panel.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size=rel(1.0))
    )
  p
  return(p)
}
