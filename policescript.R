#Packages used
library(ggplot2)
library(stm)
library(LexisNexisTools)
library(dplyr)
library(quanteda)
library(tokenizers)
library(magrittr)
library(tidytext)
library(drlib)
library(here)
#Lines 18-100 are for MODEL SELECTION which will produce
#similar topics to the ones presented in the final project
#but not exact. To load and inspect the models used, go to
#line 104 and load the .RData file and run the rest of the
#code.

set.seed(845818)

#CLEAN THE SAMPLE
#import articles
articles <- lnt_read(here("Articles"), 
                     file_type = "docx")

#remove articles with no Date
articles <- articles[-191]
articles <- articles[-240]
articles <- articles[-411]
articles <- articles[-647]
articles <- articles[-681]
articles <- articles[-780]

#custom stop words to be removed from text
mystopwords <- c("mr.", "mr", "mrs.", "ms.", "ms", "mrs", 
                 "new york times", "times reporter")

#creation of stemmed and n-grammed data-term matrix with punctuation, numbers,
#symbols, sparse words, and stop words removed
tdmtrial <- articles@articles$Article %>%
  tokens(remove_punct = TRUE, 
         remove_numbers = TRUE, 
         remove_symbols = TRUE) %>%
  tokens_remove(pattern = c(stopwords(source = "smart"), 
                            mystopwords)) %>%
  tokens_wordstem() %>%
  tokens_ngrams(n = 1:3) %>% 
  dfm() %>%
  dfm_trim(sparsity = 0.99)

#convert document-term matrix to corpus suitable for stm package 
corpustrial <- convert(tdmtrial, to = "stm")

docstrial <- corpustrial$documents
vocabtrial <- corpustrial$vocab
corpustrial$meta <- articles@meta
metatrial <- corpustrial$meta

#turning dates into numeric form for estimation of effects and plotting
corpustrial$meta$DateN <- as.numeric(corpustrial$meta$Date)
metatrial$DateN <- as.numeric(corpustrial$meta$Date)

#creating periods of time
corpustrial$meta$Periods <- case_when(as.Date(corpustrial$meta$Date) < as.Date("1991-03-03") ~ 1,
                                      as.Date(corpustrial$meta$Date) > 
                                        as.Date("1991-03-03") & as.Date(corpustrial$meta$Date)
                                      < as.Date("2014-08-10") ~ 2,
                                      as.Date(corpustrial$meta$Date) > 
                                        as.Date("2014-08-10") & as.Date(corpustrial$meta$Date)
                                      < as.Date("2020-05-26") ~ 3,
                                      as.Date(corpustrial$meta$Date) >
                                        as.Date("2020-05-25") ~ 4)

corpustrial$meta$pre <- case_when(corpustrial$meta$Periods == 1 ~ 1, 
                                 corpustrial$meta$Periods > 1 ~ 0)


corpustrial$meta$king <- case_when(corpustrial$meta$Periods == 2 ~ 1, 
                                 corpustrial$meta$Periods > 2 ~ 0,
                                 corpustrial$meta$Periods < 2 ~ 0)

corpustrial$meta$brown <- case_when(corpustrial$meta$Periods == 3 ~ 1, 
                                 corpustrial$meta$Periods > 3 ~ 0,
                                 corpustrial$meta$Periods < 3 ~ 0)

corpustrial$meta$floyd <- case_when(corpustrial$meta$Periods == 4 ~ 1, 
                                 corpustrial$meta$Periods > 4 ~ 0,
                                 corpustrial$meta$Periods < 4 ~ 0)
#inspection of top 100 words and phrases in text
topfeatures(
  tdmtrial,
  n = 100,
  decreasing = TRUE,
  scheme = c("count", "docfreq"),
  groups = NULL)

#creating diagnostic statistics for selection of topics 
kResult <- searchK(corpustrial$documents, 
                   corpustrial$vocab, 
                   K = c(5, 6, 7, 8, 9, 10, 
                         11, 12, 13, 14, 15, 
                         20, 30, 40, 50), 
                   prevalence = ~Date, 
                   data = metatrial)

#summary and plots of diagnostic statistics
kResult$results

plot(kResult$results$semcoh, 
     kResult$results$exclus,
     main = "Scatterplot of exclusivity and semantic coherence scores of structural topic models",
     xlab = "Semantic Coherence",
     ylab = "Exclusivity")

plot(kResult)

#creating multiple runs of model
runs.out <- selectModel(docstrial, 
                        vocabtrial, 
                        K = 20, 
                        prevalence = ~Date,
                        data = metatrial,
                        max.em.its = 2000,
                        seed = 845818,
                        frexw = 0.7,
                        runs = 20,
                        to.disk = T)

plotModels(runs.out)
#END OF MODEL SELECTION

#models and objects used for analysis
load(here("saved_mod.RData"))

#selected model
stm20 <- runs.out$runout[[3]]

#trial runs
stmtrial <- stm(docstrial,
                vocabtrial,
                K = 20,
                prevalence =~s(DateN) + Periods + king + brown + floyd, 
                max.em.its = 2000, 
                data = metatrial, 
                init.type = "Spectral",
                verbose = FALSE)

#plot of topic quality by exclusivity and semantic 
#coherence 
topicQuality(stmtrial, 
             docstrial)

#running multiple STM from model runs to check stability of model
stability <- multiSTM(runs.out, 
                      mass.threshold = .75,
                      reg.formula = ~Date,
                      metadata = metatrial)

#plots of model stability
plot(stability)

#plot of top 5 words in each topic
plot(stm20, n = 5)

#plot of topic associations and documents
plot.STM(stm20, 
         "hist", 
         topics = 2,
         col = "pink")

#plotting of word probabilities by topic
td_beta <- tidy(stmtrial)

options(repr.plot.width=7, 
        repr.plot.height=8, 
        repr.plot.res=100)

td_beta %>%
  group_by(topic) %>%
  top_n(10, beta) %>%
  ungroup() %>%
  mutate(topic = paste0("Topic ", topic),
         term = reorder_within(term, beta, topic)) %>%
  ggplot(aes(term, beta, fill = as.factor(topic))) +
  geom_col(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~ topic, scales = "free_y") +
  coord_flip() +
  scale_x_reordered() +
  labs(x = NULL, y = expression(beta),
       title = "Highest word probabilities for each topic",
       subtitle = "Different words are associated with different topics")

#plotting of extended word probabilities by specific topic
betaT1 <- td_beta %>%
  mutate(topic = paste0("Topic ", topic),
         term = reorder_within(term, beta, topic)) %>%
  filter(topic == "Topic 1")

betaplotT1 <- ggplot(betaT1[betaT1$beta > 0.003,], 
                     aes(term, beta, fill = as.factor(topic))) +
  geom_bar(alpha = 0.8, 
           show.legend = FALSE, 
           stat = "Identity") +
  coord_flip() +
  labs(x ="Terms", 
       y = expression(beta),
       title = "Word probabilities for Topic 17 (Police Violence)")

betaplotT1

#investigation of words by topic using different between-
#and within-topic statistics
labelTopics(stmtrial, 
            topics = c(1), 
            n = 10)

plot.STM(stm20, 
         "labels", 
         topics = c(7), 
         label = "prob", 
         n = 10, 
         width = 500)

#plotting of specific topics tool
plot(stmtrial,
     type = "labels",
     n = NULL,
     topics = 5,
     labeltype = c("prob", "frex", "lift", "score"),
     frexw = 0.5,
     main = NULL,
     xlim = NULL,
     ylim = NULL,
     xlab = NULL,
     family = "",
     width = 60,
     covarlevels = NULL,
     plabels = NULL,
     text.cex = 1,
     custom.labels = NULL,
     topic.names = NULL)

#finding exemplary articles from specific topics
plotQuote(findThoughts(stmtrial, 
                       texts = articles@articles$Headline, 
                       topics = 1, 
                       n = 1), 
          width = 45,
          text.cex = 0.75)

thoughts12 <- findThoughts(stmtrial, 
             texts = articles@meta$Headline,
             n = 12,
             topics = 12)$docs

thoughts16 <- findThoughts(stmtrial, 
                           texts = articles@meta$Headline,
                           n = 12,
                           topics = 16)$docs

par(mfrow = c(1, 2),mar = c(.5, .5, 1, .5))
plotQuote(thoughts1[c(1,4)], width = 45, main = "Topic 1")
#estimation of effect of time
stm1effect <- estimateEffect(formula =
                              1:20~s(DateN) + Periods + DateN + Periods:DateN, 
                             stmobj = stmtrial,
                             metadata = corpustrial$meta)

topicCorr(stmtrial)
plot.topicCorr(topicCorr(stmtrial))
#plot of effect of time on prevalence of specific topics

stm1effect$data$Periods %>% table()

par(bty="n",lwd=2,xaxt="n")
plot.estimateEffect(stm1effect,
                    covariate = "DateN",
                    model = stmtrial,
                    topics = stm1effect$topics[12],
                    method = "continuous",
                    xlab = "Day of the Year (January 1, 1980 to December 31, 2021)",
                    ylab = "Expected Topic Proportions",
                    main = "Coverage of Police Use of Force Incidents",
                    #moderator = "Periods",
                    #moderator.value = 1,
                    npoints = 500,
                    nsims = 500,
                    ylim = c(-0.05, 0.12),
                    xlim = c(3862, 18992),
                    linecol = "blue",
                    printlegend = F,
                    verbose.labels = T,
                    ci.level = FALSE,
                    add = F)
abline(h=0,lty=4,lwd=1,col="grey45")
abline(v=c(7731, 16292, 18408),lty=2,lwd=1,col="grey45")
par(xaxt="s")
axis(1,at=c(3652,
            7731,
            16292,
            18408,
            18992),                      
     labels=c(as.Date(3652, origin = "1970-01-01"),
              as.Date(7731, origin = "1970-01-01"),
              as.Date(16292, origin = "1970-01-01"),
              as.Date(18408, origin = "1970-01-01"), 
              as.Date(18992, origin = "1970-01-01")),
              las=1)
legend("topleft",legend=c("Charged Coverage"),
       col=c("blue"), lty=1)

plot.estimateEffect(stm1effect,
                    covariate = "DateN",
                    model = stm,
                    topics = stm1effect$topics[16],
                    method = "continuous",
                    xlab = "Day of the Year (January 1, 1980 to December 31, 2021)",
                    ylab = "Expected Topic Proportions",
                    main = "Coverage of Police Use of Force Incidents",
                    npoints = 500,
                    nsims = 500,
                    #ylim = c(-1, 1),
                    xlim = c(3862, 18992),
                    linecol = "red",
                    printlegend = F,
                    ci.level = 0,
                    add = T)

plot.estimateEffect(stm1effect,
                    covariate = "DateN",
                    model = stm,
                    topics = stm1effect$topics[12],
                    method = "continuous",
                    xlab = "Day of the Year (January 1, 1980 to December 31, 2021)",
                    ylab = "Expected Topic Proportions",
                    main = "Coverage of Police Use of Force Incidents",
                    moderator = "Periods",
                    moderator.value = 4,
                    npoints = 500,
                    nsims = 500,
                    #ylim = c(-1, 1),
                    xlim = c(3862, 18992),
                    linecol = "green",
                    printlegend = F,
                    ci.level = NULL,
                    add = T)

plot.estimateEffect(stm1effect,
                    covariate = "DateN",
                    model = stm,
                    topics = stm1effect$topics[12],
                    method = "continuous",
                    xlab = "Day of the Year (January 1, 1980 to December 31, 2021)",
                    ylab = "Expected Topic Proportions",
                    main = "Coverage of Police Use of Force Incidents",
                    moderator = "Periods",
                    moderator.value = 4,
                    npoints = 500,
                    nsims = 500,
                    #ylim = c(-1, 1),
                    xlim = c(3862, 18992),
                    linecol = "black",
                    printlegend = F,
                    ci.level = NULL,
                    add = T)



abline(h=0,lty=4,lwd=1,col="grey45")  # Put a dotted line on the y axis at 0.
abline(v=c(7731, 16292, 18408),lty=2,lwd=1,col="grey45")  # Put dotted lines 
# on the x axis at each election year.
par(xaxt="s")

#summary statistics of Date's effect on topic 17's prevalence
summary(stm1effect)
