# This is main.R
path <- "C:\\Users\\gabri\\OneDrive - Università degli Studi di Padova\\Desktop\\UNI\\Magistrale\\2022-2023\\advstat\\project\\code"
# Source the module

source(paste(path, "\\mymodule.R", sep = ""))

# Use the functions
result1 <- add_numbers(1, 2)
result2 <- multiply_numbers(3, 4)

print(result1)
print(result2)