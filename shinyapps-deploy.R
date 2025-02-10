# Edit with your shinyapps.io token information

# install rsconnect
install.packages("rsconnect")

# Set account info (Replace with actual values)
rsconnect::setAccountInfo(name='your_shinyapps_username',
                          token='your_shinyapps_token',
                          secret='your_shinyapps_secret')

# Deploy app
rsconnect::deployApp(appDir = ".", appName = "my_shiny_app")