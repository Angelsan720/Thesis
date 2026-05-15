library(tidyverse)
library(StMoMo)
library(demography)
library(zoo)                                     # Load zoo
library(glue)
library(devtools)
#install_github("PPgp/wpp2024")
library(wpp2024)
library(optimization)
normalized_distance <- function(value, upper, central, lower) {
    #case_when(
    #    value >= central & value <= upper ~ (value - central) / (upper - central),
    #    value <= central & value >= lower ~ (central - value) / (central - lower),
    #    value >= upper ~ 1,
    #    value <= lower ~ 1,
    #    TRUE ~ NA_real_
    #)
    case_when(
        value >= lower & value <= upper ~ 0,
        value > upper ~ 1,
        value < lower ~ 1,
        TRUE ~ NA_real_
    )
}
get_demo_data <- function(df, name){
    years <- df |> select(year) |> distinct() |> pull(year) |> as.numeric() |> sort() # nolint: line_length_linter.
    ages <- df |> select(age) |> distinct() |> pull(age) |> as.numeric() |> sort()
    type <- "mortality"
    label <- "PR"
    
    pop <- df |> select(c(year, age, pop)) |> 
        pivot_wider(
            names_from = age,
            values_from = pop) |>
        select(-c(year)) |> t()
    
    deaths <- df |> select(c(year, age, deaths)) |>
        pivot_wider(
            names_from = age,
            values_from = deaths) |>
        select(-c(year)) |> t()
    data <- demogdata(data = deaths/pop, pop = pop,ages = ages,years = years,type=type,label=label,name = name) # nolint
    
    
    return(data)
    
    
    
}
get_data <- function(root="data"){
    tmp_male_deaths       <- rio::import(glue('{root}/defunciones_hombres.csv')) |> select(-c("desconocido")) |>  pivot_longer( values_to = "deaths", cols = -c("year"), names_to = "age" ) |> mutate(sex="male")
    tmp_female_deaths     <- rio::import(glue('{root}/defunciones_mujeres.csv')) |> select(-c("desconocido")) |>  pivot_longer( values_to = "deaths", cols = -c("year"), names_to = "age" ) |> mutate(sex="female")
    tmp_male_population   <- rio::import(glue('{root}/poblacion_hombres.csv'))   |> pivot_longer( values_to = "pop", cols = -c("year"), names_to = "age") |>  mutate(sex="male")
    tmp_female_population <- rio::import(glue('{root}/poblacion_mujeres.csv'))   |> pivot_longer( values_to = "pop", cols = -c("year"), names_to = "age") |>   mutate(sex="female")
    df <- tmp_male_population |> bind_rows(tmp_female_population) |>
        full_join(
            tmp_male_deaths |> bind_rows( tmp_female_deaths), by=join_by(age,sex, year)) |>  filter(year>=1950)|>  filter(year<=2023) |>   mutate(year = as.numeric(year))|>   mutate(age = as.numeric(age)) 
}
compute_lt <- function(rates){
    
    #central_historic <- calc_compound_lt(rates|> filter(type=="historic") |> mutate(rate=central)) |> rename(age=x)
    central <- calc_compound_lt(rates |> mutate(rate=central)) |> rename(age=x)
    upper <- calc_compound_lt(rates |> mutate(rate=upper)) |> rename(age=x)
    lower <- calc_compound_lt(rates |> mutate(rate=lower)) |> rename(age=x)
    
    forcasted <- central |> select(year, age, ex) |>
        rename(central=ex) |>
        left_join(
            upper |> select(year, age, ex) |> rename(upper=ex),
            by=join_by(year,age)
        ) |>
        left_join(
            lower |> select(year, age, ex) |> rename(lower=ex),
            by=join_by(year,age)
        )
    
    
    return(forcasted)
}
compute_life_table3 <- function(death_rates) {
    
    m    <- death_rates$rate
    edad <- death_rates$age
    last <- length(edad)
    a <- c(rep(2.5,length(edad)))
    a[1] <- 0.045 + 2.684 * m[1] # calcula el valor de a para el primer grupo de edad
    a[2] <- 1.651 - 2.816 * m[1] # calcula el valor de a para el segundo grupo de edad
    a[last] <- NA # asigna missing al último grupo de edad
    
    n <- diff(edad)
    n[last] <- NA # asigna missing al ultimo grupo de edad
    
    q <- (n*m)/(1+(n-a)*m)
    q[last] <- 1
    p <- 1-q
    l <- 100000*c(1,cumprod(p[1:last-1]))
    l <- round(l,0)
    d <- -diff(l)
    d[last] <- l[last]
    L <- n*l - (n-a)*d # esta no es la formula que vimos enclase pero es más facil para programar
    L[last] <- l[last]/m[last] # fórmula especial para la L del último grupo de edad # nolint: line_length_linter.
    
    #i <- 1:(last-1) # `i` 'indice para accesar el proximo grupo de edad
    #L <- n*l[i+1] + a*d
    #L[last] <- l[last]/m[last]
    
    T <- rev(cumsum(rev(L)))
    e <- T/l
    a[last] <- e[last]
    
    
    life_table <- data.frame(
        x = edad,
        n = n,
        qx = q,
        lx = l,
        dx = d,
        Lx = L,
        Tx = T,
        ex = e
    )
    return (life_table)
}
calc_compound_lt<- function(df){
    return(df |> group_by(year)|> group_modify(~ compute_life_table3(.x)) |> ungroup())
}
fit_func <- function(data, link, const, wxt, oxt, head_cut=0, tail_cut=0, ages.fit = NULL, years.fit = NULL){
    #min_year <- min(data$years)
    #max_year <- max(data$years)
    #years.fit <- (min_year+tail_cut):(max_year-head_cut)
    #LCFit <- fit(lc(link = link, const = const), data = data, years.fit = years.fit, wxt = wxt, oxt = oxt,verbose = FALSE)
    LCFit <- fit(lc(link = link, const = const), data = data, wxt = wxt, oxt = oxt,verbose = FALSE, ages.fit=ages.fit, years.fit = years.fit)
    return(LCFit)
}
boot_func <- function(fit_params,nBoot, type, deathType, nSim, h,kt.method, jumpchoice, kt.lookback, oxt){
    Boot <- bootstrap(fit_params, nBoot = nBoot, type = type, deathType = deathType , oxt = oxt)
    BootSim <- simulate(Boot, nsim = nSim, h=h, jumpchoice = jumpchoice, kt.method=kt.method, kt.lookback=kt.lookback)
    
    
    tmp1 <- apply(BootSim$rates, c(1, 2), quantile, probs = 0.5) |>
        as.data.frame() |>
        rownames_to_column("age") |> pivot_longer(names_to = "year",cols = -age,values_to = "central") |>
        arrange( age,year) |> select(age,year,central)
    tmp2 <- apply(BootSim$rates, c(1, 2), quantile, probs = 0.975) |> as.data.frame() |>
        rownames_to_column("age") |>
        pivot_longer(names_to = "year",cols = -age,values_to = "lower") |>
        arrange( age,year) |> select(lower)
    tmp3 <- apply(BootSim$rates, c(1, 2), quantile, probs = 0.025) |> as.data.frame() |>
        rownames_to_column("age") |>
        pivot_longer(names_to = "year",cols = -age,values_to = "upper") |>
        arrange( age,year) |> select(upper)
    
    
    mxt <- bind_rows(
        c(tmp1,tmp2,tmp3)
    )  |> 
        mutate(age = (as.numeric(age)))|>
        mutate(year = (as.numeric(year))) |>
        arrange(age, year)
    return( mxt)
    
}

forecast_func <- function(fit_params,h,kt.method, jumpchoice, kt.lookback, oxt){
    forecast(fit_params, jumpchoice = jumpchoice, h = h, kt.method = kt.method, kt.lookback = kt.lookback, oxt = oxt,level = c(95))
}

work <- function(df, sex="female",
                 head_cut = 13,
                 tail_cut = 0,
                 rolling_window = 5,
                 const = "sum",
                 link = "log",
                 kt.method = "mrwd",
                 kt.lookback = NULL,
                 nBoot = 100,
                 type = "residual",
                 deathType = "observed",
                 nSim = 100,
                 h = 40,
                 jumpchoice = "actual",
                 wxt = NULL,
                 oxt = NULL,
                 ages.fit = NULL,
                 years.fit = NULL
                 ) {

    data <- list()
    data$Historic   <- df  |> select(c(year, age, pop, deaths)) |> mutate(rate = deaths / pop)
    data$Source   <- get_demo_data(data$Historic, name=sex) |> StMoMoData()
    data$Fit   <- fit_func(data=data$Source,   link=link, const=const, wxt, oxt, head_cut=head_cut, tail_cut=tail_cut, ages.fit=ages.fit, years.fit=years.fit)
    data$Boot   <- boot_func(data$Fit,  nBoot, type, deathType,nSim, h,kt.method, jumpchoice, kt.lookback,oxt)
    tmp = list()
    tmp$historic   <-data$Historic    |> group_by(year) |> group_modify(~compute_life_table3(.x)) |> ungroup()
    tmp$central    <-data$Boot        |> mutate(rate = central) |> group_by(year) |> group_modify(~compute_life_table3(.x)) |> ungroup()
    tmp$upper      <-data$Boot        |> mutate(rate = upper)   |> group_by(year) |> group_modify(~compute_life_table3(.x)) |> ungroup()
    tmp$lower      <-data$Boot        |> mutate(rate = lower)   |> group_by(year) |> group_modify(~compute_life_table3(.x)) |> ungroup()

    
    lt <- tmp$historic |> select(year, x, ex) |> rename(historic=ex) |>
        full_join( tmp$central |> select(year, x, ex) |> rename(central=ex), join_by(year, x)) |>
        full_join( tmp$upper |> select(year, x, ex) |> rename(upper=ex), join_by(year, x)) |>
        full_join( tmp$lower |> select(year, x, ex) |> rename(lower=ex), join_by(year, x))|>
        mutate(sex=sex)|> mutate(in_interval = ifelse(historic >= lower & historic <= upper, TRUE, FALSE)) |>
        mutate(overlap = ifelse( historic > central | historic < central  ,historic, 0)) |>
        mutate(historic_no_overlap = na_if(historic, overlap)) |>
        mutate(score = normalized_distance(historic, upper, central, lower))
    data$lt <- lt
    return(data)
    
}
work2 <- function(df, sex="female",
          head_cut = 13,
          tail_cut = 0,
          rolling_window = 5,
          const = "sum",
          link = "log",
          kt.method = "mrwd",
          kt.lookback = NULL,
          type = "residual",
          deathType = "observed",
          h = 40,
          jumpchoice = "actual",
          wxt = NULL,
          oxt = NULL,
          ages.fit = NULL,
          years.fit = NULL
) {
    
    data <- list()
    data$Historic   <- df  |> select(c(year, age, pop, deaths)) |> mutate(rate = deaths / pop)
    data$Source   <- get_demo_data(data$Historic, name=sex) |> StMoMoData()
    data$Fit   <- fit_func(data=data$Source,   link=link, const=const, wxt, oxt, head_cut=head_cut, tail_cut=tail_cut, ages.fit=ages.fit, years.fit=years.fit)
    data$Forecast   <- forecast_func(data$Fit,  h=h,kt.method=kt.method, jumpchoice=jumpchoice, kt.lookback=kt.lookback,oxt=oxt)
    return(data)
    
}
work_with_transfer <- function(df, sex="female",
                               head_cut = 13,
                               tail_cut = 0,
                               rolling_window = 5,
                               const = "sum",
                               link = "log",
                               kt.method = "mrwd",
                               kt.lookback = NULL,
                               nBoot = 100,
                               type = "residual",
                               deathType = "observed",
                               nSim = 100,
                               h = 40,
                               jumpchoice = "actual",
                               wxt = NULL,
                               oxt = NULL
) {

    data <- list()
    data$Historic         <- df  |> select(c(year, age, pop, deaths)) |> mutate(rate = deaths / pop)


    data <- list()
    data$Historic   <- df  |> select(c(year,age,pop,deaths)) |> mutate(rate = deaths/pop)

    data$Source   <- get_demo_data(data$Historic |> filter(year <= (max(year)-head_cut))|> filter(year >= (min(year)+tail_cut)), name="male"  ) |> StMoMoData()
    tail_cut = head_cut
    head_cut = 0
    
    data$Source2   <- get_demo_data(data$Historic |> filter(year <= (max(year)-head_cut))|> filter(year >= (min(year)+tail_cut)), name="male"  ) |> StMoMoData()

    
    
    data$Fit   <- fit_func(data=data$Source,   link=link, const=const, wxt, oxt, head_cut=head_cut, tail_cut=tail_cut)
    data$Fit2   <- fit_func(data=data$Source2,   link=link, const=const, wxt, oxt, head_cut=head_cut, tail_cut=tail_cut)
    
    
    data$Fit2$ax <- data$Fit$ax
    data$Fit2$bx <- data$Fit$bx
    data$Fit2$kt <- data$Fit$kt
    
    data$Boot   <- boot_func(data$Fit2,  nBoot, type, deathType,nSim, h,kt.method, jumpchoice, kt.lookback,oxt)
    tmp = list()
    tmp$historic   <-data$Historic    |> group_by(year) |> group_modify(~compute_life_table3(.x)) |> ungroup()
    tmp$central    <-data$Boot        |> mutate(rate = central) |> group_by(year) |> group_modify(~compute_life_table3(.x)) |> ungroup()
    tmp$upper      <-data$Boot        |> mutate(rate = upper)   |> group_by(year) |> group_modify(~compute_life_table3(.x)) |> ungroup()
    tmp$lower      <-data$Boot        |> mutate(rate = lower)   |> group_by(year) |> group_modify(~compute_life_table3(.x)) |> ungroup()
    
    lt <- tmp$historic |> select(year, x, ex) |> rename(historic=ex) |>
        full_join( tmp$central |> select(year, x, ex) |> rename(central=ex), join_by(year, x)) |>
        full_join( tmp$upper |> select(year, x, ex) |> rename(upper=ex), join_by(year, x)) |>
        full_join( tmp$lower |> select(year, x, ex) |> rename(lower=ex), join_by(year, x))|>
        mutate(sex=sex)|> mutate(in_interval = ifelse(historic >= lower & historic <= upper, TRUE, FALSE)) |>
        mutate(overlap = ifelse( historic > central | historic < central  ,historic, 0)) |>
        mutate(historic_no_overlap = na_if(historic, overlap)) |>
        mutate(score = normalized_distance(historic, upper, central, lower))
    data$lt <- lt
    return(data)
    }






