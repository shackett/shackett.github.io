---
title: "Building a large database of MMA fight results III: summarizing the demographics of 140,000 MMA fighters"
description: "Network-based inference of MMA fighter demographics"
category: MMA
tags: [dplyr, treemaps]
output: html_document
---

The essential units of analysis in MMA are the sport's fighters and the bouts between them. In my [previous post](https://shackett.github.io/mma/summarizingFights/), I discussed how to standardize match data so that the major categories of finishes could be easily visualized.

In this entry, I will discuss how to clean-up fighter-level data. The main goal of this post will be to determine the factors that affect fighter matchups. For some fighters, our data contains useful demographics, while for other fighters, this information is missing. One goal of this post will be to infer missing demographic data based on fighters' previous bouts. For example, if we don't know where in the world a fighter lives, we can probably determine their location with good accuracy if all their fights have been against people from the same region.



The data that we will primarily work with was collected seperately for each fighter. We will first combine fighters' data into a single table that can be more easily used.


{% highlight r %}
library(dplyr)

load("~/Desktop/MMA/software/data/raw_fight_data/fighter_output.Rdata")

kable(fighter_output[[1]]$vital_stats, row.names=F)
{% endhighlight %}



|Name            |Nickname     |Birthday   | Height| Weight|nationality |weight_class |camp             |
|:---------------|:------------|:----------|------:|------:|:-----------|:------------|:----------------|
|Andrei Arlovski |The Pit Bull |1979-02-04 | 193.04| 109.32|Belarus     |Heavyweight  |Jackson-Wink MMA |



{% highlight r %}
fighter_vital_statistics <- bind_rows(lapply(fighter_output, function(x){x$vital_stats})) %>%
  mutate(Query = names(fighter_output)) %>% tbl_df()

kable(head(fighter_vital_statistics), row.names=F)
{% endhighlight %}



|Name             |Nickname          |Birthday   | Height| Weight|nationality   |weight_class  |camp             |Query                           |
|:----------------|:-----------------|:----------|------:|------:|:-------------|:-------------|:----------------|:-------------------------------|
|Andrei Arlovski  |The Pit Bull      |1979-02-04 | 193.04| 109.32|Belarus       |Heavyweight   |Jackson-Wink MMA |/fighter/Andrei-Arlovski-270    |
|Ikuhisa Minowa   |Minowaman         |1976-01-12 | 175.26|  83.91|Japan         |Middleweight  |Kuma Gym         |/fighter/Ikuhisa-Minowa-250     |
|Kazushi Sakuraba |The Gracie Hunter |1969-07-14 | 182.88|  75.75|Japan         |Welterweight  |Laughter7        |/fighter/Kazushi-Sakuraba-84    |
|Ronda Rousey     |Rowdy             |1987-02-01 | 167.64|  61.23|United States |Bantamweight  |Team Hayastan    |/fighter/Ronda-Rousey-73073     |
|B.J. Penn        |The Prodigy       |1978-12-13 | 175.26|  65.77|United States |Featherweight |BJ Penn's MMA    |/fighter/BJ-Penn-1307           |
|Thales Trindade  |Pezao             |NA         |     NA|     NA|Brazil        |NA            |NA               |/fighter/Thales-Trindade-166551 |

### Summarizing fighter demographics: gender, weight and locale

For many fighters, our dataset contains important fighter attributes such as nationality, weight class and age. One attribute that is noticeably absent is gender. If we are interested in the factors affecting who-fights-who in MMA, then gender is definitely important.

One method that we can use to try to predict the genders of fighters is to look at the gender that is generally associated with a fighter's first name. The [gender](https://cran.rstudio.com/web/packages/gender/) provides great resources for exactly this task. Using the "ssa" method, which looks at Social Security Administration baby names, for each first name, we can determine the fraction of male babies and the fraction of female babies with the name.


{% highlight r %}
fighter_names <- fighter_vital_statistics %>% select(Name, Query) %>%
  rowwise() %>% mutate(first_name = strsplit(Name, split = " ")[[1]][1])

gender_prediction <- gender(unique(fighter_names$first_name), method = "ssa") %>% tbl_df() %>%
  select(first_name = name, Pr_male = proportion_male)

fighter_names <- fighter_names %>% left_join(gender_prediction, by = "first_name")

kable(head(fighter_names), row.names=F)
{% endhighlight %}



|Name             |Query                           |first_name | Pr_male|
|:----------------|:-------------------------------|:----------|-------:|
|Andrei Arlovski  |/fighter/Andrei-Arlovski-270    |Andrei     |  1.0000|
|Ikuhisa Minowa   |/fighter/Ikuhisa-Minowa-250     |Ikuhisa    |      NA|
|Kazushi Sakuraba |/fighter/Kazushi-Sakuraba-84    |Kazushi    |      NA|
|Ronda Rousey     |/fighter/Ronda-Rousey-73073     |Ronda      |  0.0083|
|B.J. Penn        |/fighter/BJ-Penn-1307           |B.J.       |      NA|
|Thales Trindade  |/fighter/Thales-Trindade-166551 |Thales     |  1.0000|

Classifying names works for the majority of fighters, but for 8% of fighters predicted gender is clearly ambiguous (0.1 < Pr(Male) < 0.9), while for 
11% of fighters no gender prediction is possible at all.

In addition to gender, we would like to summarize each fighter's weight class and where he/she is from. To summarize where fighters are from, nationality may be too specific for some purposes (there are 193 countries in the UN). To provide a less specific summary of fighters' locations, we can use the [countrycode](https://cran.r-project.org/web/packages/countrycode/index.html) package. Countrycode provides mappings between countries and regions. We can then also manually add the continents in which these regions fall.




{% highlight r %}
library(countrycode)

# Define regions
nationality_summary <- fighter_vital_statistics %>% count(nationality) %>%
  mutate(nation_region = countrycode(nationality, "country.name", "region"))
nationality_summary$nation_region[nationality_summary$nationality == "Taiwan"] <- "Eastern Asia"
nationality_summary$nation_region[nationality_summary$nation_region %in% c("Melanesia", "Micronesia", "Polynesia", "Australia and New Zealand")] <- "Oceania"

# Define continents
nationality_summary$nation_continent <- nationality_summary$nation_region
nationality_summary$nation_continent[nationality_summary$nation_region %in% c("Caribbean", "Northern America", "Central America")] <- "N. America"
nationality_summary$nation_continent[nationality_summary$nation_region %in% c("South America")] <- "S. America"
nationality_summary$nation_continent[grepl('Africa', nationality_summary$nation_region)] <- "Africa"
nationality_summary$nation_continent[grepl('Europe', nationality_summary$nation_region)] <- "Europe"
nationality_summary$nation_continent[grepl('Asia', nationality_summary$nation_region)] <- "Asia"

# overwrite misspelled nationalities
fighter_vital_statistics$nationality[fighter_vital_statistics$nationality %in% nationality_summary$nationality[is.na(nationality_summary$nation_region)]] <- NA

kable(head(nationality_summary), row.names = F)
{% endhighlight %}



|nationality |  n|nation_region   |nation_continent |
|:-----------|--:|:---------------|:----------------|
|165         |  1|NA              |NA               |
|185         |  1|NA              |NA               |
|Afganistan  |  1|NA              |NA               |
|Afghanistan | 98|Southern Asia   |Asia             |
|Albania     | 13|Southern Europe |Europe           |
|Algeria     | 10|Northern Africa |Africa           |

Based on our starting weight class and nationality data along with the genders that were inferred based on names, we have complete data for 49% of fighters. We can be confident in weight class (although this does change across fighters' careers). We are less confident, however, in the inferred gender because there are many cases of ambigously predicted genders.


### Network-based imputation of fighter demographics

To improve our prediction of fighters' locale, weight class and gender, we can use an additional important source of information: their matchups with other fighters. Fighters are likely to fight people who are a similar weight, the same gender and who reside in a similar location.

Our strategy is to capitalize on bout information to label fighters' unknown attributes (and gender, as this may be mispredicted by name), according to a consensus of fighters' neighbors. This approach is implemented as an iterative algorithm. During one iteration, we start with bout data (fighter-opponent) and fill in the opponent demographics. Then a consensus of opponents' demographics is taken for each fighter and the fighter's attributes are updated if they were previously unknown. After each iteration we learn more about opponents' demographics (if they were missing); this information is then passed onto fighters until we are left with only ambiguous attributes.



{% highlight r %}
# load bout data
all_bouts <- readRDS("~/Desktop/MMA/software/data/processed_fight_data/MMA_all_bouts.Rds")

fighter_vital_network_V <- fighter_vital_statistics %>% select(Query, weight_class, nationality) %>%
  left_join(fighter_names, by = "Query") %>%
  left_join(nationality_summary %>% select(-n), by = "nationality") %>%
  mutate(Pr_male_impute = Pr_male,
         weight_class_impute = weight_class,
         nationality_impute = nationality,
         nation_region_impute = nation_region,
         nation_continent_impute = nation_continent)

kable(head(fighter_vital_network_V), row.names = F)
{% endhighlight %}



|Query                           |weight_class  |nationality   |Name             |first_name | Pr_male|nation_region    |nation_continent | Pr_male_impute|weight_class_impute |nationality_impute |nation_region_impute |nation_continent_impute |
|:-------------------------------|:-------------|:-------------|:----------------|:----------|-------:|:----------------|:----------------|--------------:|:-------------------|:------------------|:--------------------|:-----------------------|
|/fighter/Andrei-Arlovski-270    |Heavyweight   |Belarus       |Andrei Arlovski  |Andrei     |  1.0000|Eastern Europe   |Europe           |         1.0000|Heavyweight         |Belarus            |Eastern Europe       |Europe                  |
|/fighter/Ikuhisa-Minowa-250     |Middleweight  |Japan         |Ikuhisa Minowa   |Ikuhisa    |      NA|Eastern Asia     |Asia             |             NA|Middleweight        |Japan              |Eastern Asia         |Asia                    |
|/fighter/Kazushi-Sakuraba-84    |Welterweight  |Japan         |Kazushi Sakuraba |Kazushi    |      NA|Eastern Asia     |Asia             |             NA|Welterweight        |Japan              |Eastern Asia         |Asia                    |
|/fighter/Ronda-Rousey-73073     |Bantamweight  |United States |Ronda Rousey     |Ronda      |  0.0083|Northern America |N. America       |         0.0083|Bantamweight        |United States      |Northern America     |N. America              |
|/fighter/BJ-Penn-1307           |Featherweight |United States |B.J. Penn        |B.J.       |      NA|Northern America |N. America       |             NA|Featherweight       |United States      |Northern America     |N. America              |
|/fighter/Thales-Trindade-166551 |NA            |Brazil        |Thales Trindade  |Thales     |  1.0000|South America    |S. America       |         1.0000|NA                  |Brazil             |South America        |S. America              |



{% highlight r %}
# Look at all opponents
fighter_vital_network_E <- all_bouts %>% select(Fighter_link, Opponent_link) %>% unique()
fighter_vital_network_E <- rbind(fighter_vital_network_E, fighter_vital_network_E %>% select(Opponent_link = Fighter_link, Fighter_link = Opponent_link)) %>% unique()

# Add known attributes of query fighter (fighter link)
fighter_vital_network_E <- fighter_vital_network_E %>% left_join(fighter_vital_network_V %>% select(Fighter_link = Query, Pr_male, weight_class, nationality, nation_region, nation_continent), by = "Fighter_link")
# Add consensus estimate of opponent attributes (opponent_link)
fighter_vital_network_E <- fighter_vital_network_E %>% left_join(fighter_vital_network_V %>% select(Opponent_link = Query, Pr_male_impute:nation_continent_impute), by = "Opponent_link")

kable(head(fighter_vital_network_E), row.names = F)
{% endhighlight %}



|Fighter_link                 |Opponent_link                      | Pr_male|weight_class |nationality |nation_region  |nation_continent | Pr_male_impute|weight_class_impute |nationality_impute |nation_region_impute |nation_continent_impute |
|:----------------------------|:----------------------------------|-------:|:------------|:-----------|:--------------|:----------------|--------------:|:-------------------|:------------------|:--------------------|:-----------------------|
|/fighter/Andrei-Arlovski-270 |/fighter/Stipe-Miocic-39537        |       1|Heavyweight  |Belarus     |Eastern Europe |Europe           |             NA|Heavyweight         |United States      |Northern America     |N. America              |
|/fighter/Andrei-Arlovski-270 |/fighter/Frank-Mir-2329            |       1|Heavyweight  |Belarus     |Eastern Europe |Europe           |         0.9953|Heavyweight         |United States      |Northern America     |N. America              |
|/fighter/Andrei-Arlovski-270 |/fighter/Travis-Browne-16785       |       1|Heavyweight  |Belarus     |Eastern Europe |Europe           |         0.9930|Heavyweight         |United States      |Northern America     |N. America              |
|/fighter/Andrei-Arlovski-270 |/fighter/Antonio-Silva-12354       |       1|Heavyweight  |Belarus     |Eastern Europe |Europe           |         0.9914|Heavyweight         |Brazil             |South America        |S. America              |
|/fighter/Andrei-Arlovski-270 |/fighter/Brendan-Schaub-33926      |       1|Heavyweight  |Belarus     |Eastern Europe |Europe           |         0.9960|Heavyweight         |United States      |Northern America     |N. America              |
|/fighter/Andrei-Arlovski-270 |/fighter/Andreas-Kraniotakes-30848 |       1|Heavyweight  |Belarus     |Eastern Europe |Europe           |         0.9920|Heavyweight         |Germany            |Western Europe       |Europe                  |



{% highlight r %}
# determine the top match from neighbors
get_max <- function(x){
  if(all(is.na(x))){
    return(c(NA_character_)) 
  }else{
    x %>% table() %>% sort(decreasing = T) %>% names() %>% first() %>% return()
  }
}

# so that we can track missing values as the algorithm proceeds
count_na <- function(x){
 sum(is.na(x)) 
}

# track number of missing values throughout iteration

iter_count <- 0
track_fill_in <- fighter_vital_network_V %>% ungroup() %>% summarize_each(funs(count_na), Pr_male_impute:nation_continent_impute) %>% mutate(iteration = iter_count)

repeat{
  iter_count <- iter_count+1
  
  # infer demographics based on neighborhood
  # for weight class and locations, find the most frequent neighbor class
  # for gender [0,1] average Pr(male)
  
  update_fighter_vitals <- fighter_vital_network_E %>%
    group_by(Fighter_link) %>%
    summarize(Pr_male = Pr_male[1],
              weight_class = weight_class[1],
              nationality = nationality[1],
              nation_region = nation_region[1],
              nation_continent = nation_continent[1],
              Pr_male_impute = ifelse(all(is.na(Pr_male_impute)), NA, mean(Pr_male_impute, na.rm = T)),
              weight_class_impute = get_max(weight_class_impute),
              nationality_impute = get_max(nationality_impute),
              nation_region_impute = get_max(nation_region_impute),
              nation_continent_impute = get_max(nation_continent_impute))
  
  track_fill_in <- rbind(track_fill_in,
                         update_fighter_vitals %>% summarize_each(funs(count_na), Pr_male_impute:nation_continent_impute) %>% mutate(iteration = iter_count)
  )
  
  if(iter_count > 9){
    break
  }else{
    
    # overwrite known demographics during each iteration
    update_fighter_vitals <- update_fighter_vitals %>%
      mutate(Pr_male_impute = ifelse(is.na(Pr_male), Pr_male_impute, Pr_male),
             weight_class_impute = ifelse(is.na(weight_class), weight_class_impute, weight_class),
             nationality_impute = ifelse(is.na(nationality), nationality_impute, nationality),
             nation_region_impute = ifelse(is.na(nation_region), nation_region_impute, nation_region),
             nation_continent_impute = ifelse(is.na(nation_continent), nation_continent_impute, nation_continent))
    
    fighter_vital_network_E <- fighter_vital_network_E %>% select(Fighter_link:nation_continent) %>%
      left_join(update_fighter_vitals %>% select(Opponent_link = Fighter_link, Pr_male_impute:nation_continent_impute),
                by = "Opponent_link")
    
  }
}

fighter_vital_consensus <- update_fighter_vitals %>% rowwise() %>%
  mutate(Gender = if(is.na(Pr_male_impute)){NA_character_}else if(Pr_male_impute < 0.2){"Female"}else if(Pr_male_impute > 0.8){"Male"}else{NA_character_},
         Weight_class = ifelse(is.na(weight_class), weight_class_impute, weight_class),
         Nationality = ifelse(is.na(nationality), nationality_impute, nationality),
         Nation_region = ifelse(is.na(nation_region), nation_region_impute, nation_region),
         Nation_continent = ifelse(is.na(nation_continent), nation_continent_impute, nation_continent)) %>%
  select(Query = Fighter_link, Gender:Nation_continent)
  
fighter_vital_statistics <- fighter_vital_statistics %>% select(Query, Name:Weight, Camp = camp) %>% left_join(fighter_vital_consensus, by = "Query")
{% endhighlight %}

As the imputation algorithm proceeds, we are able to fill in a large amount of missing information. This method should accurately predict unknown genders, weight classes and the general locale of fighters (inferred nationalities are more uncertain in many cases).



{% highlight r %}
kable(track_fill_in)
{% endhighlight %}



| Pr_male_impute| weight_class_impute| nationality_impute| nation_region_impute| nation_continent_impute| iteration|
|--------------:|-------------------:|------------------:|--------------------:|-----------------------:|---------:|
|          15890|               55810|              42317|                42317|                   42317|         0|
|           8719|               26693|              20429|                20429|                   20429|         1|
|            820|               13393|               9595|                 9595|                    9595|         2|
|            688|               11144|               7775|                 7775|                    7775|         3|
|            671|               10680|               7412|                 7412|                    7412|         4|
|            670|               10576|               7338|                 7338|                    7338|         5|
|            670|               10551|               7319|                 7319|                    7319|         6|
|            670|               10546|               7312|                 7312|                    7312|         7|
|            670|               10544|               7308|                 7308|                    7308|         8|
|            670|               10542|               7307|                 7307|                    7307|         9|
|            670|               10542|               7307|                 7307|                    7307|        10|

### Summary of fighter demographics

As a quick summary of cleaned-up fighter-specific data, we can look at summaries of two major types of fighter demographics: fighter weight classes and where fighters live.

First we will look at weight classes. Because weight classes are ordered, with Atomweight the lightest and Super Heavyweight the heaviest, weight categories were sorted according to their logical progression and the frequency of fighters in each weight class was determined


{% highlight r %}
library(ggplot2)

weight_summary <- fighter_vital_statistics %>%
  count(Weight_class) %>%
  filter(!is.na(Weight_class)) %>%
  mutate(Weight_class = ordered(Weight_class,
                                levels = c("Atomweight", "Strawweight", "Flyweight", "Bantamweight",
                                           "Featherweight", "Lightweight", "Welterweight", "Middleweight",
                                           "Light Heavyweight", "Heavyweight", "Super Heavyweight"))) %>%
  arrange(Weight_class)

# ggplot2 plotting theme
hist_theme <- theme(axis.title = element_text(color = "black", size = 25),
                    panel.background = element_rect(fill = "gray93"),
                    panel.border = element_rect(fill = NA, color = "black"),
                    panel.grid.minor = element_blank(),
                    panel.grid.major.x = element_blank(),
                    panel.grid.major.y = element_line(color = "gray70", size = 1),
                    axis.text = element_text(size = 18, hjust = 0.5, vjust = 0.5),
                    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                    axis.ticks = element_line(size = 1),
                    axis.ticks.length = unit(0.15, "cm"))

ggplot(weight_summary, aes(x = Weight_class, y = n)) +
  geom_bar(stat = "identity", fill = "firebrick") +
  hist_theme +
  scale_y_continuous("Counts") + scale_x_discrete("Weight Class")
{% endhighlight %}

<img src="/figure/source/2016-05-13-summarizingFighters/unnamed-chunk-3-1.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

Looking at weight classes: most fighters compete at intermediate weight classes; the most common weight classes are Lightweight and Welterweight.

The geographical summaries that we generated: nationality, region and continent are a nested hierarchy with nations residing within regions and regions within continents. To make use of this hierarchy for visualization we need to make sure that this nested structure is maintained for all fighters such that nations can only belong to one region and regions can only belong to one continent.


{% highlight r %}
locale_summary <- fighter_vital_statistics %>%
  count(Nation_continent, Nation_region, Nationality)

locale_summary <- locale_summary %>%
  filter(!is.na(Nation_continent) & !is.na(Nation_region) & !is.na(Nationality))

# assigning each nation to unique region and region to a continent
continent_summary <- locale_summary %>%
  count(Nation_continent, Nation_region, wt = n) %>%
  group_by(Nation_region) %>%
  summarize(Continent = Nation_continent[which.max(n)][1])

region_summary <- locale_summary %>%
  count(Nation_region, Nationality, wt = n) %>%
  group_by(Nationality) %>%
  summarize(Region = Nation_region[which.max(n)][1])

locale_summary <- locale_summary %>%
  left_join(region_summary, by = "Nationality") %>%
  count(Nation_continent, Region, Nationality, wt = n) %>%
  left_join(continent_summary, by = c("Region" = "Nation_region")) %>%
  count(Continent, Region, Nationality, wt = n) %>%
  ungroup

kable(head(locale_summary %>% arrange(desc(n)), 8), row.names = F)
{% endhighlight %}



|Continent  |Region           |Nationality    |     n|
|:----------|:----------------|:--------------|-----:|
|N. America |Northern America |United States  | 57170|
|S. America |South America    |Brazil         | 18553|
|Europe     |Northern Europe  |United Kingdom |  8190|
|Asia       |Eastern Asia     |Japan          |  5928|
|N. America |Northern America |Canada         |  4618|
|Europe     |Eastern Europe   |Russia         |  4478|
|N. America |Northern America |USA            |  2436|
|Europe     |Eastern Europe   |Poland         |  2361|

Having established a hierarchy of fighter locations, summarized based on the number of fighters residing in each country, nested within geographical regions and continents, we can visualize this hierarchy using a treemap. Essentially we will turn each country into a rectangle with area proportional to the number of fighter's from that country. These rectangles will be nested according to the continent-region-nationality hierarchy.

To allow for clearer discrimination of regions (which will be especially useful for later analyses) we can shade similar countries with similar colors. Conventionally, we would only color one level of our hierarchy by either choosing unique colors for continents or regions or nations. If we only colored continents, we would only have a very coarse summary of geograophy. If we colored nations then similar country might receive very different shades of color.

Using a hierarchy to guide the choice of colors seems natural when such information is present. Here, I chose a general color for each continent, shades of that color for regions, and more specific colors for each country.


{% highlight r %}
library(treemap)

source("~/Desktop/MMA/manuscripts/blogs/Coloring_hierarchical_data/hierarchical_color_lib.R")

available_colors <- extract_color_space(hmax = 355, lmin = 30, cmin = 30)

country_colors <- identify_color_hierarchy(locale_summary, available_colors, weight_column = "n") %>%
  filter(Tier == "Nationality") %>% select(Nationality = Category, Color)

locale_summary <- locale_summary %>% ungroup %>% left_join(country_colors, by = "Nationality")

treemap(locale_summary, 
index = c("Continent", "Region", "Nationality"),
vSize = "n",
vColor = "Color", type = "color",
title = "Locations of MMA fighters",
fontsize.title = 25,
fontsize.labels = c(0, 0, 15))
{% endhighlight %}

<img src="/figure/source/2016-05-13-summarizingFighters/unnamed-chunk-5-1.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" style="display: block; margin: auto;" />

By coloring and organizing fighter according to their nationality we can see that Mixed Martial Arts is truely a world-wide phenomenon, although its popularity is greatest in the USA, Europe, Brazil and Japan.

In my next post, I will demonstrate how the wide geographical distribution of MMA fighters and other factors such as weight class and gender structure the matchups between fighters. This structure emerges directly from the matchups between fighters and can be readily visualized by treating MMA fights as an undirected graph.
