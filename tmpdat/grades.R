d1 = read.csv("/Users/Chris/Dropbox/Mac/Downloads/grades.csv") %>%
  mutate(Email = email)


d2 = read.csv("/Users/Chris/Dropbox/Mac/Downloads/POL 201 SP25 002_GradesExport_2025-03-17-15-52.csv")

# Join the two
# dataframes on the "Student ID" column
library(dplyr)
d3 = left_join(d1, d2, by = "Email")

write.csv(d3, "/Users/Chris/Dropbox/Mac/Downloads/grades_combined.csv", row.names = FALSE)
