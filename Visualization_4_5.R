##### Visualization for Assumption (4) and Assumption (5) #####

source("library.R")

power_table_3<-get(load("data/power_table_3.Rdata"))

### POWER TABLE ###

Power_3 <-power_table_3 %>%
  mutate(n = factor(n)) %>%
  pivot_longer(c("reject_mv_mv", "reject_mv_sn", "reject_mv_st", "reject_sn_st"), names_to = "test", values_to = "reject") %>%
  group_by(mu, sigma, test, n) %>%
  summarize(power = mean(reject))



### VISUALIZATION ###
ui <- fluidPage(
  theme = bslib::bs_theme(
    bg = "grey", 
    fg = "white", 
    base_font = "Source Sans Pro"),
  fg = "white",
  titlePanel("Comparing Mean Vectors from Two Populations"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Power study comapring he mean vectors from two populations with different covariance matrices different and without the normality assumption."),
      
      radioButtons("Test", "Choose a test",
                   #levels(as.factor(Power_3$test)),
                   choiceNames = c("normal vs normal", " normal vs. skew-normal","normal vs. skew-t", "skew-normal vs. skew-t"),
                   choiceValues = levels(as.factor(Power_3$test))),
      
      checkboxGroupInput("N", "Choose a n", 
                         choices = as.numeric(levels(as.factor(Power_3$n))),
                         selected = as.numeric(levels(as.factor(Power_3$n)))),
      
      sliderTextInput(
        inputId = "Sigma",
        label = "Choose k such that Sigma1 = k Sigma2",
        choices = as.numeric(levels(as.factor(Power_3$sigma))),
        selected = range(Power_3$sigma)[1],
        animate = TRUE,
        grid =TRUE
      )
    ),
    mainPanel(plotOutput("out"))
  )
)

server <- function(input, output) {
  output$out <- 
    renderPlot({ 
      filter(Power_3, n %in% input$N, sigma == input$Sigma, test == input$Test)  %>%
        ggplot() + geom_line(aes(x = mu, y = power, by = test, linetype = n), color = "blue") +
        geom_hline(yintercept = 0.01) +
        labs(title= "Power study of two Mean Vectors")+
        ylab("Estimated Power") +
        xlab(expression(mu))
    })
}

shinyApp(ui, server)



## TYPE-I ERROR ###
knitr::kable(arrange(Power_3, n) %>% filter(mu ==0))


arrange(Power_3, n) %>% filter(mu ==0) %>%
  group_by(n) %>%
  summarize(mean_power = mean(power))


arrange(Power_3, n) %>% filter(mu ==0) %>%
  group_by(test) %>%
  summarize(mean_power = mean(power))

arrange(Power_3, n) %>% filter(mu ==0) %>%
  group_by(sigma) %>%
  summarize(mean_power = mean(power))
