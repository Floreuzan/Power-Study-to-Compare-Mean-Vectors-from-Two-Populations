##### Visualization for Assumption (5) #####

source("library.R")

power_table_2<-get(load("data/power_table_2.Rdata"))

### POWER TABLE ###

Power_2 <-power_table_2 %>%
  mutate(n = factor(n)) %>%
  group_by(mu, sigma, n) %>%
  summarize(power = mean(reject))

### VISUALIZATION ###

ui <- fluidPage(theme = bslib::bs_theme(
  bg = "grey", 
  fg = "white", 
  base_font = "Source Sans Pro"),
  titlePanel("Comparing Mean Vectors from Two Populations"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Power study comapring he mean vectors from two populations with different covariance matrices."),
      
      checkboxGroupInput("N", "Choose a n", 
                         choices = as.numeric(levels(as.factor(Power_2$n))),
                         selected = as.numeric(levels(as.factor(Power_2$n)))),
      sliderTextInput(
        inputId = "Sigma",
        label = "Choose k such that Sigma1 = k Sigma2",
        choices = as.numeric(levels(as.factor(Power_2$sigma))),
        selected = range(Power_2$sigma)[1],
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
      filter(Power_2, n %in% input$N, sigma == input$Sigma)  %>%
        ggplot() + geom_line(aes(x = mu, y = power, linetype = n), color = "blue") +
        geom_hline(yintercept = 0.01) +
        labs(title= "Power study of two Mean Vectors")+
        ylab("Estimated Power") +
        xlab(expression(mu))
    })
}
shinyApp(ui, server)


### TYPE-I ERROR ###
knitr::kable(arrange(Power_2, n) %>% filter(mu ==0))

arrange(Power_2, n) %>% filter(mu ==0) %>%
  group_by(n) %>%
  summarize(mean_power = mean(power))

arrange(Power_2, n) %>% filter(mu ==0) %>%
  group_by(sigma) %>%
  summarize(mean_power = mean(power))
