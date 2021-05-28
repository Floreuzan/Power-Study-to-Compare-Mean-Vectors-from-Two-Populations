##### Visualization for Assumption (4) #####

source("library.R")

power_table_1<-get(load("data/power_table_1.Rdata"))

### POWER TABLE ###
Power_1 <-power_table_1 %>%
  mutate(n = factor(n)) %>%
  pivot_longer(c("reject_mv_mv", "reject_sn_sn","reject_st_st", "reject_mv_sn", "reject_mv_st", "reject_sn_st"), names_to = "test", values_to = "reject") %>%
  group_by(mu, test, n) %>%

### VISUALIZATION ###
ui <- fluidPage(theme = bslib::bs_theme(
  bg = "grey", 
  fg = "white", 
  base_font = "Source Sans Pro"),
  titlePanel("Comparing Mean Vectors from Two Populations"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Power study comparing he mean vectors from two populations without the normality assumption.."),
      
      radioButtons("Test", "Choose a test",
                   #levels(as.factor(Power_1$test)),
                   choiceNames = c("normal vs normal","skew-normal vs skew-normal", "skew-t vs skew-t", "normal vs. skew-normal","normal vs. skew-t ", "skew-normal vs skew-t"),
                   choiceValues = levels(as.factor(Power_1$test))),
      
      checkboxGroupInput("N", "Choose a n", 
                         choices = as.numeric(levels(as.factor(Power_1$n))),
                         selected = as.numeric(levels(as.factor(Power_1$n))))
    ),
    mainPanel(plotOutput("out"))
  )
)

server <- function(input, output) {
  output$out <- 
    renderPlot({ 
      filter(Power_1, n %in% input$N, test == input$Test)  %>%
        ggplot() + geom_line(aes(x = mu, y = power, by = test, linetype = n), color = "blue") +
        geom_hline(yintercept = 0.01) +
        labs(title= "Power study of two Mean Vectors")+
        ylab("Estimated Power") +
        xlab(expression(mu))
    })
}
shinyApp(ui, server)


### TYPE-I ERROR ###
knitr::kable(arrange(Power_1, n) %>% filter(mu ==0))

arrange(Power_1, n) %>% filter(mu ==0) %>%
  group_by(n) %>%
  summarize(mean_power = mean(power))

arrange(Power_1, n) %>% filter(mu ==0) %>%
  group_by(test) %>%
  summarize(mean_power = mean(power))