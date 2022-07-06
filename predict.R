pred <- dat %>%
  mutate(p_hat = predict(fit, newdata = dat, type = "response"))

# logit prediction

pred <- pred %>%
  mutate(logit = log(p_hat / (1 - p_hat))) %>%
  gather(key = "predictors", value = "predictor.value", -logit)

pred$p_hat[is.na(pred$p_hat)] <- 0

#plots all predictor values

ggplot(pred, aes(logit, predictor.value)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess") + 
  facet_wrap(~predictors, scales = "free_y")

#plots individual predictor value

pred %>%
  filter(predictors=='k7_AntigensSelfTest') %>%
  ggplot(aes(logit, predictor.value)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess")
