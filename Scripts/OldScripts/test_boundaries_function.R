#test boundaries 
Arguments


#inner	
#logical. If TRUE, "inner" boundaries are returned, else "outer" boundaries are returned

#classes	
#character. Logical. If TRUE all different values are (after rounding) distinguished, as well as NA. If FALSE (the default) only edges between NA and non-NA cells are considered

#directions	
#integer. Which cells are considered adjacent? Should be 8 (Queen's case) or 4 (Rook's case)

#falseval	
#numeric. The value to use for cells that are not a boundary and not NA


r <- rast(nrows=18, ncols=36, xmin=0)
plot(r)
r[150:250] <- 1
r[251:450] <- 2
plot(r)
forest_mask <- r == 2  # This creates a logical (TRUE/FALSE) mask where TRUE is for forest (2)
plot(forest_mask)


#now find cells where cells border a different class values in at least one NSEW direction 
class_difference <- boundaries(r,inner= TRUE, directions = 4, classes=TRUE)
plot(class_difference)

#now find only class difference cells that are forest 
x <- class_difference*forest_mask
plot(x)
plot(r)
