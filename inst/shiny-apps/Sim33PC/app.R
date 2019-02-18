#
# This is a Shiny web application intended to demonstrate the concepts of:
#   Norris DC. Precautionary Coherence Unravels Dose Escalation Designs.
#     bioRxiv. December 2017:240846. doi:10.1101/240846
#     https://www.biorxiv.org/content/early/2017/12/29/240846
#
# TODO:
# 1. Arrange the inputs in a row beneath the D3 app:
# (a) Median & CV of MTDi
# (b) [Simulate] button
# (c) # doses; max & min dose
# (d) # periods
# 2. Activate the [Simulation] button
# 3. Feed doses (min, max, #) to de.sim()
# 4. Set # periods to plot
# 5. Set # participants to enroll
# 6. Should I hand-craft an HTML myself, with DIV container, etc?
# 7. Do dynamic updating of D3 axes, according to selected parameters?
# ≈ Note that while some inputs (like dose level selection) would seem
#   to wipe the slate clean, others allow for more graded updates. For
#   example, changing the number of periods and/or enrollment could
#   update the existing OXDS plot. It may even help clarify the overall
#   visualization if I try to decouple these various kinds of update.
# 8. Highlight relevant portions of D3 viz when inputs get focus
# 9. Modify Gamma cumulative distribution plot dynamically
# ≈ Note that, if CV is limited to a fixed number of levels, the D3
#   viz might even 'predraw' generic versions of all possible dists,
#   stretching and shifting these as needed (if SVG supports this).

library(shiny)
library(shinyWidgets)
library(r2d3)
library(DTAT)

# signif(10^seq(0, 1, 1/40), 3):
forty <- c(1.00, 1.06, 1.12, 1.19, 1.26, 1.33, 1.41, 1.50, 1.58, 1.68,
           1.78, 1.88, 2.00, 2.11, 2.24, 2.37, 2.51, 2.66, 2.82, 2.99,
           3.16, 3.35, 3.55, 3.76, 3.98, 4.22, 4.47, 4.73, 5.01, 5.31,
           5.62, 5.96, 6.31, 6.68, 7.08, 7.50, 7.94, 8.41, 8.91, 9.44)
# Manual nudges to 'painfully close' dose levels:
forty[c(4,5,7,15,17,20)] <- c(1.2, 1.25, 1.4, 2.25, 2.5, 3.0)
forty[c(21,24:29)] <- c(3.15, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0)
forty[c(32:35,37:39)] <- c(6.0, 6.3, 6.7, 7.0, 8.0, 8.4, 8.9)
# Create 4 decades
decadence <- c(forty/10, forty, forty*10, forty*100, 1000)

ui <- fluidPage(
   
   d3Output("d3"),
   
   hr(),

   fluidRow(
     column(4,
            sliderInput("D"
                         ,"Number of doses"
                         ,min = 3
                         ,max = 9
                         ,value = 7)
            ,sliderInput("N"
                         ,"Number of participants"
                         ,min = 18
                         ,max = 60
                         ,step = 3
                         ,value = 24)
     ),
     column(4,
            shinyWidgets::sliderTextInput("Drange","Dose range",
                  choices=decadence,
                  selected=c(0.25, 1.88), grid = F)
            ,prettyRadioButtons("Dunit"
                         ,"Dose Units"
                         ,c("ng"="ng",
                            "µg"="µg",
                            "mg"="mg",
                            "g"="g")
                         ,selected = "mg"
                         ,inline = TRUE)
            ,prettyRadioButtons("Dper"
                         ,NULL
                         ,c("/m²"="/m²",
                            "/kg"="/kg",
                            "absolute"="")
                         ,selected = "/kg"
                         ,inline = TRUE)
     ),
     column(4,
            shinyWidgets::sliderTextInput("MTDi50","Median MTDi",
                                          choices=decadence,
                                          selected=c(1.0), grid = F)
            ,sliderInput("CV"
                         ,"CV(MTDi)"
                         ,min = 20
                         ,max = 150
                         ,step = 5
                         ,value = 70
                         ,post = "%")
            ,actionBttn("simButton"
                        ,"Simulate Your 3+3/PC Trial")
     )
   ),
   hr(),
   tags$div(
      class = "header", checked = NA,
      tags$p(tags$b("See:"),
             "Norris DC. Precautionary Coherence Unravels Dose Escalation Designs",
             tags$i("bioRxiv"), "December 2017:240846.",
             tags$a(href="https://www.biorxiv.org/content/early/2017/12/29/240846",
                    "doi:10.1101/240846.")
             )
   )
)

server <- function(input, output) {
   
   output$d3 <- renderD3({
      CV <- as.numeric(input$CV)/100
      shape <- CV^-2
      scale <- as.numeric(input$MTDi50)/qgamma(0.5, shape=shape)
      start.dose <- input$Drange[1]
      dose.jump <- (input$Drange[2]/input$Drange[1]) ^ (1/(input$D-1)) - 1
      N <- input$N
      input$simButton # 'take a dependency' on the simButton
      # Let's check what values get passed in:
      cat(paste("start.dose =", start.dose, "\n"))
      cat(paste("highest.dose =", input$Drange[2], "\n"))
      cat(paste("dunit =", paste(input$Dunit, input$Dper, sep=''), "\n"))
      cat(paste("dose.jump =", dose.jump, "\n"))
      cat(paste("n.doses =", input$D, "\n"))
      cat(paste("CV =", CV, "\n"))
      cat(paste("mean_mtd =", shape*scale, "\n"))
      cat(paste("median MTD =", input$MTDi50, "\n"))
      doses <- start.dose*(1+dose.jump)^(0:(input$D-1))
      trial <- new("DE", doses=doses,
                   MTDi=rgamma(N, shape=shape, scale=scale),
                   units="mg/kg")
      trial <- titration(trial, periods=10)
      plot(trial)
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

