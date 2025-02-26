#' Make a posterior predictive check plot
#'
#' @param dt Output of `make_posterior_predictive_checks()`
#' @param outdir Directory to save generated plots
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' dt.po <- extract_posterior_predictions(fit, dt)
#' dt.po <- make_posterior_predictive_check(dt.po)
#'
#' plot_posterior_predictive_checks(dt.po)
#' }
plot_posterior_predictive_checks <- function(dt, outdir=NA){
  plt.ppd_check <- ggplot(dt) +
    geom_tile(aes(x = age, y = alter_age, fill = factor(inside.CI))) +
    labs(x = "Age of contacting individuals",
         y = "Age of contacted individuals",
         fill = "Inside 95% CI" ) +
    coord_equal() +
    facet_grid(alter_gender_label ~ gender_label) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(color=NA, fill = "transparent")
    )

  if(!is.na(outdir)){
    ggsave(file.path(outdir, "figures", "ppd_check.png"), plot = plt.ppd_check)
  } else {
    return(plt.ppd_check)
  }
}

#' Make a predicted contacts plot
#'
#' @param dt Output of `extract_posterior_predictions()`
#' @param outdir Directory to save generated plots
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' dt.po <- extract_posterior_predictions(fit, dt)
#' plot_predicted_contacts(dt.po)
#' }
plot_predicted_contacts <- function(dt, outdir=NA){
  plt <- ggplot(dt) +
    geom_tile(aes(x = age, y = alter_age, fill = M)) +
    labs(x = "Age of contacting individuals", y = "Age of contacted individuals", fill = "Predicted contacts" ) +
    coord_equal() +
    facet_grid(alter_gender_label ~ gender_label) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    guides(fill = guide_colorbar(title.position = "top", barwidth = 10)) +
    scale_fill_continuous(type = "viridis", option="B", na.value = "white",
                          limits = c(min(dt$y), max(dt$y))) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(color=NA, fill = "transparent")
    )

  if(!is.na(outdir)){
    ggsave(file.path(outdir, "figures", "predicted_contacts.png"), plot = plt)
  } else {
    return(plt)
  }
}

#' Make a predicted contact intensities plot
#'
#' @param dt Output of `extract_posterior_predictions()`
#' @param outdir Directory to save generated plots
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' dt.po <- extract_posterior_predictions(fit, dt)
#' plot_posterior_intensities(dt.po)
#' }
plot_posterior_intensities <- function(dt, outdir=NA, new_hh=FALSE){
  if (new_hh){
    p <- ggplot(dt) +
      geom_tile(aes(x = age, y = alter_age, fill = M)) +
      labs(x = "Participants' age", y = "Contacts' age", fill = "Contact rate" ) +
      coord_equal() +
      facet_grid(paste(alter_gender, "(Contacts)") ~ paste(gender, "(Participants)")) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      viridis::scale_fill_viridis(na.value="white", option="H") +
      theme_bw() +
      theme(
        legend.position = "bottom",
        strip.background = element_rect(color=NA, fill = "transparent"),
        text = element_text(size = 8)) 
      # scale_fill_continuous(labels = scales::label_number(scale = 10000, suffix = "k")) 
    
    if(!is.na(outdir)){
      ggsave(file.path(outdir, "figures", "rate_matrices.png"), plot = p)
    } else {
      return(p)
    }
  }
  else{
    p <- ggplot(dt) +
    geom_tile(aes(x = age, y = alter_age, fill = intensity_M)) +
    labs(x = "Participants' age", y = "Contacts' age", fill = "Contact intensity" ) +
    coord_equal() +
    facet_grid(paste(alter_gender, "(Contacts)") ~ paste(gender, "(Participants)")) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    viridis::scale_fill_viridis(na.value="white", option="H") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(color=NA, fill = "transparent")
    )

  if(!is.na(outdir)){
    ggsave(file.path(outdir, "figures", "intensity_matrices.png"), plot = p)
  } else {
    return(p)
  }
  }
}

#' Make a predicted contact rates plot
#'
#' @param dt Output of `extract_posterior_predictions()`
#' @param dt.strata Output of `extract_posterior_predictions()` (age stratified)
#' @param outdir Directory to save generated plots
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' dt.po <- extract_posterior_predictions(fit, dt)
#' plot_predicted_rates(dt.po)
#' }
plot_predicted_rates <- function(dt, outdir=NA){
  plt <- ggplot(dt) +
    geom_tile(aes(x = age, y = alter_age, fill = cntct_rate_predict)) +
    labs(x = "Age of contacting individuals", y = "Age of contacted individuals", fill = "Predicted contact rates" ) +
    coord_equal() +
    facet_grid(alter_gender_label ~ gender_label) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    guides(fill = guide_colorbar(title.position = "top", barwidth = 10)) +
    scale_fill_continuous(type = "viridis", na.value = "white",
                          limits = c(min(dt$cntct_rate), max(dt$cntct_rate))) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(color=NA, fill = "transparent")
    )

  if(!is.na(outdir)){
    ggsave(file.path(outdir, "figures", "predicted_rates.png"), plot = plt)
  } else {
    return(plt)
  }
}

plot_alpha <- function(dt, outdir=NA){
  p <- ggplot(dt) +
    geom_tile(aes(x = age, y = alter_age, fill = alpha_agg)) +
    labs(x = "Participants' age", y = "Contacts' age", fill = "Alpha" ) +
    facet_grid( paste(alter_gender, "(Contacts)") ~ paste(gender, "(Participants)") ) +
    coord_equal() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    viridis::scale_fill_viridis(na.value="white", option="H") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(color=NA, fill = "transparent")
    )
  
  if(!is.na(outdir)){
    ggsave(file.path(outdir, "alpha_contact_intensities.png"), plot = p)
  }
  return(p)

}

#' Plot the predicted contact rates sliced at different contacting ages
#'
#' @param dt data.table with extracted posterior predictions
#' @param age.cut a vector specifying the ages to slice at
#' @param .gender gender of contacting individual
#' @param .alter_gender gender of contacted individual
#' @param outdir directory to save plot
#'
#' @return
#' @export
#'
#' @examples
plot_sliced_intensities <- function(dt, age.cut = c(20, 40, 60), outdir=NA, new_hh=FALSE, new_hh_intensity=FALSE)
{
  if (new_hh){
    tmp <- dt[age %in% age.cut]
    tmp[, age_label := paste("Participants' age:", age)]
    tmp[, comb := paste(gender, "to", alter_gender)]
    
    p <- ggplot(tmp, aes(x=alter_age)) +
      geom_stepribbon(aes(ymin = CL, ymax = CU), alpha=0.3) +
      geom_step(aes(y = M)) +
      labs(x="Contacts' age", y="Alpha") +
      facet_grid(age_label ~ comb) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        strip.background = element_rect(color=NA, fill = "transparent"),
        strip.text = element_text(face = "bold")
      )
    if(!is.na(outdir)){
      ggsave(file.path(outdir, "figures", "sliced_alphas.png"), plot = p, width=7, height=6)
    }
    
    return(p)
  }
  
  else if (new_hh_intensity){
      tmp <- dt[age %in% age.cut]
      tmp[, age_label := paste("Participants' age:", age)]
      tmp[, comb := paste(gender, "to", alter_gender)]
      
      p <- ggplot(tmp, aes(x=alter_age)) +
        geom_stepribbon(aes(ymin = CL, ymax = CU), alpha=0.3) +
        geom_step(aes(y = M)) +
        labs(x="Contacts' age", y="Contact intensity") +
        facet_grid(age_label ~ comb) +
        theme_bw() +
        theme(
          legend.position = "bottom",
          strip.background = element_rect(color=NA, fill = "transparent"),
          strip.text = element_text(face = "bold")
        )
      if(!is.na(outdir)){
        ggsave(file.path(outdir, "figures", "sliced_intensities.png"), plot = p, width=7, height=6)
      }
      
      return(p)
    }
  else{
    tmp <- dt[age %in% age.cut]
    tmp[, age_label := paste("Participants' age:", age)]
    tmp[, comb := paste(gender, "to", alter_gender)]
    
    p <- ggplot(tmp, aes(x=alter_age)) +
      geom_stepribbon(aes(ymin = intensity_CL, ymax = intensity_CU), alpha=0.3) +
      geom_step(aes(y = intensity_M)) +
      labs(x="Contacts' age", y="Contact intensity") +
      facet_grid(age_label ~ comb) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        strip.background = element_rect(color=NA, fill = "transparent"),
        strip.text = element_text(face = "bold")
      )
    if(!is.na(outdir)){
      ggsave(file.path(outdir, "figures", "sliced_intensities.png"), plot = p, width=7, height=6)
    }
    
    return(p)
  }
}

plot_marginal_intensities <- function(dt, outdir=NA, new_hh=FALSE, rate=FALSE){
  
  if (new_hh){
    if (rate){
      
      p <- ggplot(dt, aes(age, rate_M)) +
        geom_stepribbon(aes(ymin = rate_CL, ymax = rate_CU, fill = comb), alpha=0.3) +
        geom_step(aes(color = comb)) +
        scale_y_continuous(limits = c(0, NA)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_color_brewer(palette = "Set1") +
        scale_fill_brewer(palette = "Set1") +
        labs(y="Marginal contact rate", x="Age of contacting individual",
             color="Gender combination", fill="Gender combination") +
        theme_bw() +
        theme(
          legend.position = "right",
          plot.background = element_rect(fill='transparent', color=NA),
          strip.background = element_rect(color=NA, fill = "transparent"),
          legend.background = element_rect(fill="transparent", color=NA)
        )
      
    }
    else{
      p <- ggplot(dt, aes(age, alpha_M)) +
        geom_stepribbon(aes(ymin = alpha_CL, ymax = alpha_CU, fill = comb), alpha=0.3) +
        geom_step(aes(color = comb)) +
        scale_y_continuous(limits = c(0, NA)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_color_brewer(palette = "Set1") +
        scale_fill_brewer(palette = "Set1") +
        labs(y="Marginal contact intensity", x="Age of contacting individual",
             color="Gender combination", fill="Gender combination") +
        theme_bw() +
        theme(
          legend.position = "right",
          plot.background = element_rect(fill='transparent', color=NA),
          strip.background = element_rect(color=NA, fill = "transparent"),
          legend.background = element_rect(fill="transparent", color=NA)
        )
    }
    
  }else{
    p <- ggplot(dt, aes(age, intensity_M)) +
      geom_stepribbon(aes(ymin = intensity_CL, ymax = intensity_CU, fill = comb), alpha=0.3) +
      geom_step(aes(color = comb)) +
      scale_y_continuous(limits = c(0, NA)) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_color_brewer(palette = "Set1") +
      scale_fill_brewer(palette = "Set1") +
      labs(y="Marginal contact intensity", x="Age of contacting individual",
           color="Gender combination", fill="Gender combination") +
      theme_bw() +
      theme(
        legend.position = "right",
        plot.background = element_rect(fill='transparent', color=NA),
        strip.background = element_rect(color=NA, fill = "transparent"),
        legend.background = element_rect(fill="transparent", color=NA)
      )
  }

  if(!is.na(outdir)){
    if (new_hh){
      if (rate){
        ggsave(file.path(outdir, "figures", "marginal_rates.png"), plot = p, height = 3, width = 7)
      }else{
        ggsave(file.path(outdir, "figures", "marginal_alphas.png"), plot = p, height = 3, width = 7)
      }
    }else{
      ggsave(file.path(outdir, "figures", "marginal_intensities.png"), plot = p, height = 3, width = 7)
    }
  }

  return(p)
}

#' Plots the empirical contact intensity patterns
#'
#' @param dt A data.table with the aggregated counts (ideally with completed zeros)
#' @param outdir Directory path the save the outputs
#'
#' @return ggplot
#' @export
#'
#' @examples
plot_empirical_intensities <- function(dt, outdir=NA){
  dt[, intensity := y/part]
  dt[is.nan(intensity), intensity := 0]

  p <- ggplot(dt, aes(age, alter_age_strata)) +
    geom_tile(aes(fill = intensity)) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    viridis::scale_fill_viridis(na.value = "white", option="H") +
    labs(x="Participants' age", y="Contacts' age strata", fill="Contact intensity") +
    facet_grid( paste(alter_gender, "(Contacts)") ~ paste(gender, "(Participants)") ) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      legend.position = "bottom", 
      text = element_text(size = 5),
      strip.background = element_rect(color=NA, fill = "transparent")
    )

  if(!is.na(outdir)){
    ggsave(file.path(outdir, "figures", "empirical_intensities.png"), plot = p)
  }

  return(p)
}

