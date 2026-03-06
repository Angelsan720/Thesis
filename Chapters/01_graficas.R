library(tidyverse)
library(StMoMo)
library(demography)
library(zoo)
library(glue)
library(ggfan)
fig_width <- 8
fig_height <- 4
get_demo_data <- function(df){
    years <- df |> select(year) |> distinct() |> pull(year) |> as.numeric() |> sort()
    ages <- df |> select(age) |> distinct() |> pull(age) |> as.numeric() |> sort()
    type <- "mortality"
    label <- "PR"
    name <- "female"
    pop <- df |> select(c(year, age, pop, sex))|> filter(sex==name) |> 
        pivot_wider(
            names_from = age,
            values_from = pop) |>
        select(-c(year,sex)) |> t()
    
    deaths <- df |> select(c(year, age, deaths, sex))|> filter(sex==name) |>
        pivot_wider(
            names_from = age,
            values_from = deaths) |>
        select(-c(year,sex)) |> t()
    
    female_data <- demogdata(data = deaths/pop, pop = pop,ages = ages,years = years,type=type,label=label,name = name)
    
    ###################################
    name <- "male"
    pop <- df |> select(c(year, age, pop, sex))|> filter(sex==name) |> 
        pivot_wider(
            names_from = age,
            values_from = pop) |>
        select(-c(year,sex)) |> t()
    
    deaths <- df |> select(c(year, age, deaths, sex))|> filter(sex==name) |>
        pivot_wider(
            names_from = age,
            values_from = deaths) |>
        select(-c(year,sex)) |> t()
    
    male_data <- demogdata(data = deaths/pop, pop = pop,ages = ages,years = years,type=type,label=label,name = name)
    
    
    return(list(male=male_data, female=female_data))
    
    
    
}
get_data <- function(){
    tmp_male_deaths       <- rio::import('data/defunciones_hombres.csv') |> select(-c("desconocido")) |>  pivot_longer( values_to = "deaths", cols = -c("year"), names_to = "age" ) |> mutate(sex="male")
    tmp_female_deaths     <- rio::import('data/defunciones_mujeres.csv') |> select(-c("desconocido")) |>  pivot_longer( values_to = "deaths", cols = -c("year"), names_to = "age" ) |> mutate(sex="female")
    tmp_male_population   <- rio::import('data/poblacion_hombres.csv')   |> pivot_longer( values_to = "pop", cols = -c("year"), names_to = "age") |>  mutate(sex="male")
    tmp_female_population <- rio::import('data/poblacion_mujeres.csv')   |> pivot_longer( values_to = "pop", cols = -c("year"), names_to = "age") |>   mutate(sex="female")
    df <- tmp_male_population |> bind_rows(tmp_female_population) |>
        full_join(
            tmp_male_deaths |> bind_rows( tmp_female_deaths), by=join_by(age,sex, year)) |>  filter(year>=1950)|>  filter(year<=2021) |>   mutate(year = as.numeric(year))|>   mutate(age = as.numeric(age)) 
}
linear_interpolate_vectorized <- function(sequence, i, j) {
    # Input validation
    if (i <= 0 || j >= length(sequence) || i >= j) {
        stop("Invalid indices: must satisfy 0 < i < j < n")
    }
    
    n <- length(sequence)
    
    # Extract the values at indices i and j
    y_i <- sequence[i]
    y_j <- sequence[j]
    
    # Create indices for interpolation
    indices <- (i + 1):(j - 1)
    
    # Calculate interpolated values using vectorized operations
    sequence[indices] <- y_i + ((y_j - y_i) / (j - i)) * (indices - i)
    
    return(sequence)
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
    L[last] <- l[last]/m[last] # fórmula especial para la L del último grupo de edad
    
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
func <- function(
        data,
        link, const,
        nBoot, type, deathType,
        nsim, h,kt.method, jumpchoice, kt.lookback){
    
    # wxt <- genWeightMat(data$ages,  data$years, zeroCohorts = c((max(data$years)-16):max(data$years)))
    # oxt <- rowMeans(log(data$Dxt / data$Ext))
    # LCFit <- fit(lc(link = link, const = const), data = data, wxt = wxt,oxt = oxt)
    # LCFit <- fit(lc(link = link, const = const), data = data, oxt = oxt)
    
    LCFit <- fit(lc(link = link, const = const), data = data)
    # LCFit <- fit(lc(link = link, const = const), data = data, wxt = wxt)
    
    #LCFitFor <- forecast(LCFit, h = h, jumpchoice = "actual")
    LCFitBoot <- bootstrap(LCFit, nBoot = nBoot, type = type, deathType = deathType )
    # LCFitBootSim <- simulate(LCFitBoot, nsim = nsim, h=h, jumpchoice = jumpchoice, kt.method=kt.method, oxt=oxt, kt.lookback=kt.lookback)
    LCFitBootSim <- simulate(LCFitBoot, nsim = nsim, h=h, jumpchoice = jumpchoice, kt.method=kt.method, kt.lookback=kt.lookback)
    #LCFitSim <- simulate(LCFit, nsim = n_sim, h = h, jumpchoice = "actual")
    
    
    tmp1 <- (LCFit$Dxt / LCFit$Ext) |> as.data.frame() |> rownames_to_column("age")|> pivot_longer(names_to = "year",cols = -age,values_to = "central") |> arrange( age,year)
    #tmp2 <- LCFitFor$rates |>as.data.frame() |> rownames_to_column("age")|> pivot_longer(names_to = "year",cols = -age,values_to = "central") |> arrange( age,year) |> select(age,year,central)
    tmp2 <- apply(LCFitBootSim$rates, c(1, 2), quantile, probs = 0.5) |> as.data.frame() |>
        rownames_to_column("age") |> pivot_longer(names_to = "year",cols = -age,values_to = "central") |>
        arrange( age,year) |> select(age,year,central)
    tmp3 <- apply(LCFitBootSim$rates, c(1, 2), quantile, probs = 0.975) |> as.data.frame() |>
        rownames_to_column("age") |> pivot_longer(names_to = "year",cols = -age,values_to = "lower") |>
        arrange( age,year) |> select(lower)
    tmp4 <- apply(LCFitBootSim$rates, c(1, 2), quantile, probs = 0.025) |> as.data.frame() |>
        rownames_to_column("age") |> pivot_longer(names_to = "year",cols = -age,values_to = "upper") |>
        arrange( age,year) |> select(upper)
    
    
    mxt <- bind_rows(
        c(tmp2,tmp3,tmp4)
    ) |>
        mutate(type="forecast") |>
        bind_rows(tmp1|>  mutate(type="historic"))|>
        mutate(age = (as.numeric(age)))|>
        mutate(year = (as.numeric(year))) |> arrange(age, year)
    structure(
        list(
            mxt = mxt,
            LCFit = LCFit,
            LCFitBoot = LCFitBoot,
            LCFitBootSim = LCFitBootSim
        )
    )
    
}
plot_lt <- function(df, title,age_filter){
    
    df|>filter(age==age_filter) |> ggplot(aes(y=central, x=year , colour = type, group = type)) +
        geom_ribbon(aes(ymax = upper, ymin= lower), fill = "grey70")+
        geom_line() +
        labs(
            title=title,
            x = "Year",
            y = "Life Expectancy"
        ) +
        theme_bw() +
        theme(legend.position="none") 
    
}
plot_mx <- function(df, title, age_filter){
    
    df|>filter(age==age_filter) |> ggplot(aes(y=log(central), x=year , colour = type, group = type)) +
        geom_ribbon(aes(ymax = log(upper), ymin= log(lower)), fill = "grey70")+
        geom_line() +
        labs(
            title=title,
            x = "Year",
            y = "Log Death Rate"
        ) +
        theme_bw() +
        theme(legend.position="none") 
    
}
calc_compound_lt<- function(df){
    return(df |> group_by(year)|> group_modify(~ compute_life_table3(.x)) |> ungroup())
}
compute_lt <- function(rates){
    
    central_historic <- calc_compound_lt(rates|> filter(type=="historic") |> mutate(rate=central)) |> rename(age=x)
    central <- calc_compound_lt(rates|> filter(type=="forecast") |> mutate(rate=central)) |> rename(age=x)
    upper <- calc_compound_lt(rates|> filter(type=="forecast") |> mutate(rate=upper)) |> rename(age=x)
    lower <- calc_compound_lt(rates|> filter(type=="forecast") |> mutate(rate=lower)) |> rename(age=x)
    
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
    lt <- forcasted |> mutate(type="forecast") |> bind_rows(
        central_historic |> select(year, age, ex) |> rename(central=ex) |> mutate(type="historic")
    ) |> arrange(age, year)
    
    return(lt)
}
plot_residuals_extremes <- function(data, title){
    years <- data$years
    ages <- data$ages
    df <- data.frame(
        age = rep(ages, length(years)),
        year = rep(years, each=length(ages)),
        residual = as.vector(data$residuals)
    )
    ggplot(df, aes(x=year, y=as.factor(age), fill=abs(residual)>2)) +
        geom_tile() +
        scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
        labs(title=title, x="Year", y="Age", fill="Residual") +
        theme_bw()
}
plot_residuals <- function(data_male, data_female){
    years <- data_male$years
    ages <- data_male$ages
    df <- data.frame(
        age = rep(ages, length(years)),
        year = rep(years, each=length(ages)),
        residual = as.vector(data_male$residuals)
    ) |> mutate(sex="male") |> bind_rows(
        data.frame(
            age = rep(ages, length(years)),
            year = rep(years, each=length(ages)),
            residual = as.vector(data_female$residuals)
        ) |> mutate(
            sex="female"
        ))
    
    ggplot(df, aes(x=year, y=as.factor(age), fill=residual)) +
        geom_tile() +
        scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
        labs(title="Residuals fit by sex: 1950-2022", x="Year", y="Age", fill="Residual") +
        theme_bw()+
        facet_wrap(~sex)
}
plot_params <- function(ax_male, bx_male, kt_male,ax_female, bx_female, kt_female){
    #kt_frame <- data.frame(
    #    year = mx_male$LCFit$years,
    #    kt = kt
    #)
    df_male <- data.frame(age = mx_male$LCFit$ages) |>
        mutate(ax  = ax_male) |>
        mutate(bx  = bx_male) |>
        mutate(sex = "male")
    df_female <- data.frame(age = mx_female$LCFit$ages) |>
        mutate(ax  = ax_female) |>
        mutate(bx  = bx_female) |>
        mutate(sex = "female")
    frame <- df_male |> bind_rows(df_female)
    
    
    plt1<- ggplot(frame, aes(x=age,y=ax,color=sex, group=sex)) +
        geom_line(aes(), size=1) +
        labs(title="LC parameter ax", x="Age", y="Ax") +
        theme_bw()+
        scale_color_manual(values = c("female" = "red", "male" = "blue"))
    plt2<- ggplot(frame, aes(x=age,y=bx,color=sex, group=sex)) +
        geom_line(aes(), size=1) +
        labs(title="LC parameter bx", x="Age", y="Bx") +
        theme_bw()+
        scale_color_manual(values = c("female" = "red", "male" = "blue"))
    kt_frame <- bind_rows(
        kt_male   |> as.data.frame() |> mutate(sex="male")   |> filter(sex=="male"),
        kt_female |> as.data.frame() |> mutate(sex="female") |> filter(sex=="female")
    )|> pivot_longer(names_to = "year", cols = -sex) |> rename(kt=value) |> mutate(year=as.numeric(year))
    
    plt3<- ggplot(kt_frame, aes(x=year,y=kt,color=sex, group=sex)) +
        geom_line(aes(), size=1) +
        labs(title="LC parameter kt", x="Year", y="Kt") +
        theme_bw()+
        scale_color_manual(values = c("female" = "red", "male" = "blue"))
    return(structure(  list(  ax_plot = plt1,  bx_plot = plt2,  kt_plot = plt3 )  ))
    
}


data <- list()


tmp <- get_data() |> group_by(age, sex) |>
    arrange(year) |> 
    ungroup() |> 
    filter(!is.na(deaths))
const <- "sum"
link <- "log"
kt.method <-"mrwd" 
kt.lookback = NULL
nBoot <- 50
type="residual"
deathType="observed"
nsim <- 50
h<- 40
jumpchoice="actual"
data$Source <- get_demo_data(tmp)
n100 <- 100
n1000 <- 1000
mx_female            <- func(data=data$Source$female |> StMoMoData(), link=link, const=const, nBoot=nBoot,    type=type, deathType=deathType, nsim=nsim, h=h, jumpchoice=jumpchoice, kt.method=kt.method, kt.lookback=kt.lookback)
mx_male              <- func(data=data$Source$male   |> StMoMoData(), link=link, const=const, nBoot=nBoot,    type=type, deathType=deathType, nsim=nsim, h=h, jumpchoice=jumpchoice, kt.method=kt.method, kt.lookback=kt.lookback)
#mx_boot_00100_female <- func(data=data$Source$female |> StMoMoData(), link=link, const=const, nBoot=n100,     type=type, deathType=deathType, nsim=n100, h=h, jumpchoice=jumpchoice, kt.method=kt.method, kt.lookback=kt.lookback)
#mx_boot_00100_male   <- func(data=data$Source$male   |> StMoMoData(), link=link, const=const, nBoot=n100,     type=type, deathType=deathType, nsim=n100, h=h, jumpchoice=jumpchoice, kt.method=kt.method, kt.lookback=kt.lookback)
#mx_boot_01000_female <- func(data=data$Source$female |> StMoMoData(), link=link, const=const, nBoot=n1000,    type=type, deathType=deathType, nsim=n1000, h=h, jumpchoice=jumpchoice, kt.method=kt.method, kt.lookback=kt.lookback)
#mx_boot_01000_male   <- func(data=data$Source$male   |> StMoMoData(), link=link, const=const, nBoot=n1000,    type=type, deathType=deathType, nsim=n1000, h=h, jumpchoice=jumpchoice, kt.method=kt.method, kt.lookback=kt.lookback)
#mx_boot_10000_female <- func(data=data$Source$female |> StMoMoData(), link=link, const=const, nBoot=10000,   type=type, deathType=deathType, nsim=nsim, h=h, jumpchoice=jumpchoice, kt.method=kt.method, kt.lookback=kt.lookback)
#mx_boot_10000_male   <- func(data=data$Source$male   |> StMoMoData(), link=link, const=const, nBoot=10000,   type=type, deathType=deathType, nsim=nsim, h=h, jumpchoice=jumpchoice, kt.method=kt.method, kt.lookback=kt.lookback)


residuals_male <- residuals(mx_male$LCFit)
residuals_female <- residuals(mx_female$LCFit)

ax_male <- mx_male$LCFit$ax
bx_male <- mx_male$LCFit$bx
kt_male <- mx_male$LCFit$kt
ax_female <- mx_female$LCFit$ax
bx_female <- mx_female$LCFit$bx
kt_female <- mx_female$LCFit$kt

plt <- plot_residuals(residuals_male, residuals_female)
ggsave("graphs/01_residuals.png",plot = plt, width=fig_width, height=fig_height)

plts <- plot_params(ax_male, bx_male, kt_male,ax_female, bx_female, kt_female)
ggsave("graphs/01_ax.png",plot = plts$ax_plot, width=fig_width, height=fig_height)
ggsave("graphs/01_bx.png",plot = plts$bx_plot, width=fig_width, height=fig_height)
ggsave("graphs/01_kt.png",plot = plts$kt_plot, width=fig_width, height=fig_height)

boot_male   <- mx_male$LCFitBoot
boot_female <- mx_female$LCFitBoot
plot_boot_params <- function(male_boot,female_boot, nBoots=nBoot){
    
    tmp_male <- boot_male$model$ax |> as.data.frame() 
    colnames(tmp_male) <- "ax_observed" 
    tmp_male<-tmp_male |> mutate(sex="male") |> filter(sex=="male") |> rownames_to_column("age") |> mutate(age=as.numeric(age))
    
    tmp_female <- boot_female$model$ax |> as.data.frame()
    colnames(tmp_female) <- "ax_observed" 
    tmp_female<-tmp_female |> mutate(sex="female") |> filter(sex=="female") |> rownames_to_column("age") |> mutate(age=as.numeric(age))
    
    ax_observed <- bind_rows(
        tmp_male,
        tmp_female
    )
    
    
    tmp_male <- boot_male$model$bx |> as.data.frame() 
    colnames(tmp_male) <- "bx_observed" 
    tmp_male<-tmp_male |> mutate(sex="male") |> filter(sex=="male") |> rownames_to_column("age") |> mutate(age=as.numeric(age))
    
    tmp_female <- boot_female$model$bx |> as.data.frame() 
    colnames(tmp_female) <- "bx_observed" 
    tmp_female<-tmp_female |> mutate(sex="female") |> filter(sex=="female") |> rownames_to_column("age") |> mutate(age=as.numeric(age))
    
    bx_observed <- bind_rows(
        tmp_male,
        tmp_female
    )

    kt_observed <- bind_rows(
        boot_male$model$kt   |> as.data.frame() |> mutate(sex="male")   |> filter(sex=="male"),
        boot_female$model$kt |> as.data.frame() |> mutate(sex="female") |> filter(sex=="female")
    )|> pivot_longer(names_to = "year", cols = -sex) |> rename(kt_real=value) |> mutate(year=as.numeric(year))
        

    
    
    
    
    ax_male   <- data.frame(age=mx_female$LCFit$ages)
    ax_female <- data.frame(age=mx_female$LCFit$ages)
    bx_male   <- data.frame(age=mx_female$LCFit$ages)
    bx_female <- data.frame(age=mx_female$LCFit$ages)
    kt        <- data.frame()
    for(i in 1:nBoots){
        #
        tmp <- as.data.frame(male_boot$bootParameters[[i]][["ax"]])
        colnames(tmp) <- "ax"
        tmp <- tmp |> rownames_to_column("age") |> mutate(boot=i) |> mutate(age=as.numeric(age))
        ax_male <- bind_rows(ax_male, tmp)
        #
        tmp <- as.data.frame(male_boot$bootParameters[[i]][["bx"]])
        colnames(tmp) <- "bx"
        tmp <- tmp |> rownames_to_column("age") |> mutate(boot=i) |> mutate(age=as.numeric(age))
        bx_male <- bind_rows(bx_male, tmp)
        #
        tmp <- as.data.frame(female_boot$bootParameters[[i]][["ax"]])
        colnames(tmp) <- "ax"
        tmp <- tmp |> rownames_to_column("age") |> mutate(boot=i) |> mutate(age=as.numeric(age))
        ax_female <- bind_rows(ax_female, tmp)
        #
        tmp <- as.data.frame(female_boot$bootParameters[[i]][["bx"]])
        colnames(tmp) <- "bx"
        tmp <- tmp |> rownames_to_column("age") |> mutate(boot=i) |> mutate(age=as.numeric(age))
        bx_female <- bind_rows(bx_female, tmp)
        #
        kt_frame <- bind_rows(
            male_boot$bootParameters[[i]][["kt"]]   |> as.data.frame() |> mutate(sex="male")   |> filter(sex=="male"),
            female_boot$bootParameters[[i]][["kt"]] |> as.data.frame() |> mutate(sex="female") |> filter(sex=="female")
        )|> pivot_longer(names_to = "year", cols = -sex) |> rename(kt=value) |> mutate(year=as.numeric(year)) |> mutate(boot=i)
        kt <- bind_rows(kt, kt_frame)        
    }
    kt <- kt|> left_join(kt_observed, by=join_by(year, sex)) |> mutate(delta_kt = kt - kt_real)
    ax_male <- ax_male|> mutate(sex="male")
    ax_female <- ax_female|> mutate(sex="female")
    ax<- ax_male |> bind_rows(ax_female)
    
    
    bx_male <- bx_male|> mutate(sex="male")
    bx_female <- bx_female|> mutate(sex="female")
    bx<- bx_male |> bind_rows(bx_female)
    
    df <- ax |> left_join(bx, by=join_by(age, boot, sex))
    df <- df |> left_join(ax_observed, by=join_by(age, sex))|> left_join(bx_observed, by=join_by(age, sex)) |> filter(!is.na(bx)) |> mutate(delta_ax = ax - ax_observed) |> mutate(delta_bx = bx - bx_observed)
    
    ###################### fan plot
    #plt1 <- ggplot(df, aes(x=age, y=ax)) +
    #    geom_fan() +
    #    theme_bw() +
    #    facet_wrap(~sex, ncol=1)+
    #    labs(title="Bootstrap ax", x="Age", y="Ax")
    #plt2 <- ggplot(df, aes(x=age, y=bx)) +
    #    geom_fan() +
    #    theme_bw() +
    #    facet_wrap(~sex, ncol=1)+
    #    labs(title="Bootstrap bx", x="Age", y="Bx")
    #plt3 <- ggplot(kt, aes(x=year, y=kt)) +
    #    geom_fan() +
    #    theme_bw() +
    #    facet_wrap(~sex, ncol=1)+
    #    labs(title="Bootstrap kt", x="year", y="kt")
    ###################### scatter plot
    plt4 <- ggplot(df, aes(x=as.factor(age), y=delta_ax)) +
        geom_point() +
        theme_bw() +
        facet_wrap(~sex, ncol=1)+
        labs(title="Bootstrap residuals for ax", x="Age", y="Ax")
    plt5 <- ggplot(df, aes(x=as.factor(age), y=delta_bx)) +
        geom_point() +
        theme_bw() +
        facet_wrap(~sex, ncol=1)+
        labs(title="Bootstrap residuals for bx", x="Age", y="Bx")
    plt6 <- ggplot(kt, aes(x=year, y=delta_kt)) +
        geom_point() +
        theme_bw() +
        facet_wrap(~sex, ncol=1)+
        labs(title="Bootstrap residuals for kt", x="year", y="kt")
    ###################### delta scatter plot
    ax_quantiles <- df |> group_by(age, sex) |> summarise(
        q_0025 = quantile(delta_ax, probs=0.025),
        q_0500 = quantile(delta_ax, probs=0.5),
        q_0975 = quantile(delta_ax, probs=0.975)
    )
    bx_quantiles <- df |> group_by(age, sex) |> summarise(
        q_0025 = quantile(delta_bx, probs=0.025),
        q_0500 = quantile(delta_bx, probs=0.5),
        q_0975 = quantile(delta_bx, probs=0.975)
    )
    kt_quantiles <- kt |> group_by(year, sex) |> summarise(
        q_0025 = quantile(delta_kt, probs=0.025),
        q_0500 = quantile(delta_kt, probs=0.5),
        q_0975 = quantile(delta_kt, probs=0.975)
    )
    
    return(structure(list(boot_ax_plot=plt4, boot_bx_plot=plt5, boot_kt_plot=plt6, quantiles = structure(list(ax=ax_quantiles, bx=bx_quantiles, kt=kt_quantiles)))))
}
get_quantiles <- function(sequence, quantiles=c(0.025,0.5, 0.975)){
    quantile_values <- quantile(sequence, probs = quantiles)
    return(quantile_values)
}
plts<-plot_boot_params(boot_male, boot_female)

ggsave("graphs/01_boot_ax.png",plot = plts$boot_ax_plot, width=fig_width, height=fig_height)
ggsave("graphs/01_boot_bx.png",plot = plts$boot_bx_plot, width=fig_width, height=fig_height)
ggsave("graphs/01_boot_kt.png",plot = plts$boot_kt_plot, width=fig_width, height=fig_height)

rio::export(plts$quantiles$kt, "graphs/01_boot_residuals_kt.csv","csv")
rio::export(plts$quantiles$bx, "graphs/01_boot_residuals_bx.csv","csv")
rio::export(plts$quantiles$ax, "graphs/01_boot_residuals_ax.csv","csv")

get_sim_kt <- function(male_sim, female_sim){
    sim <- male_sim
    kts   <- sim$kt.s$sim
    q_0975_male <- apply(kts, c(1, 2), quantile, probs = 0.975) |> as.vector()
    q_0500_male <- apply(kts, c(1, 2), quantile, probs = 0.050) |> as.vector()
    q_0025_male <- apply(kts, c(1, 2), quantile, probs = 0.025) |> as.vector()
    
    sim <- female_sim
    kts   <- sim$kt.s$sim
    q_0975_female <- apply(kts, c(1, 2), quantile, probs = 0.975) |> as.vector()
    q_0500_female <- apply(kts, c(1, 2), quantile, probs = 0.050) |> as.vector()
    q_0025_female <- apply(kts, c(1, 2), quantile, probs = 0.025) |> as.vector()
    
    
    years <- sim$kt.s$years
    df<- bind_rows(
        data.frame(
            year = years,
            q_0975 = q_0975_male,
            q_0500 = q_0500_male,
            q_0025 = q_0025_male,
            sex="male"
        ),
        data.frame(
            year = years,
            q_0975 = q_0975_female,
            q_0500 = q_0500_female,
            q_0025 = q_0025_female,
            sex="female"
        )
    )
}
sim_kt <- get_sim_kt(mx_male$LCFitBootSim, mx_female$LCFitBootSim)
rio::export(sim_kt, "graphs/01_simulated_kt.csv","csv")

plt_sim_kt <- sim_kt |> ggplot(aes(x=year)) +
    geom_line(aes(y=q_0975, color = sex), linetype = "dashed") +
    geom_line(aes(y=q_0500, color = sex)) +
    geom_line(aes(y=q_0025, color = sex), linetype = "dashed") +
    scale_color_manual(values = c("male" = "blue", "female"="red"))+
    theme_bw() +
    labs(title="Simulated kt by sex 2022-2072", x="Year", y="Kt")
ggsave("graphs/01_simulated_kt.png",plot = plt_sim_kt, width=fig_width, height=fig_height)
#rio::export(df, "graphs/01_life expted_kt.csv","csv")
df <- bind_rows(mx_male$mx |> mutate(sex="male"),mx_female$mx |> mutate(sex="female"))

rio::export(df, "graphs/01_deathrates.csv","csv")

compute_lt2 <- function(rates){
    
    central_historic <- calc_compound_lt(rates|> filter(type=="historic") |> mutate(rate=central)) |> rename(age=x) |> mutate(type = "historic")
    central <- calc_compound_lt(rates|> filter(type=="forecast") |> mutate(rate=central)) |> rename(age=x)|> mutate(type = "central")
    upper <- calc_compound_lt(rates|> filter(type=="forecast") |> mutate(rate=upper)) |> rename(age=x)|> mutate(type = "upper")
    lower <- calc_compound_lt(rates|> filter(type=="forecast") |> mutate(rate=lower)) |> rename(age=x)|> mutate(type = "lower")
    
    
    
    return(bind_rows(central_historic,central,upper,lower))
}

df<- df |> group_by(sex) |> group_modify(~ compute_lt2(.x))
rio::export(df, "graphs/01_lifetables.csv","csv")







