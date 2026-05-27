km_res <- plot_km_target_source_combined(
  out_file = NULL,
  ncol = 3
)

p_km <- km_res$plot +
  labs(title = "Kaplan-Meier curves: target and selected sources") +
  guides(
    color = guide_legend(ncol = 8, byrow = TRUE),
    linetype = guide_legend(ncol = 8, byrow = TRUE)
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = grid::unit(1.2, "cm"),
    legend.text = element_text(size = 7),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 8),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 10)
  )

p_cumh <- panel_plot +
  plot_annotation(title = "Baseline Cumulative Hazard Inference: Cox-SieveTL vs CoxPH Target") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 8),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 10)
  )

combined_two_panel_vertical <- cowplot::plot_grid(
  p_km,
  p_cumh,
  labels = c("A", "B"),
  label_size = 18,
  label_fontface = "bold",
  ncol = 1,
  rel_heights = c(1, 1)
)

ggsave(
  file.path(plot_dir, "KM_and_cumh_two_panel_vertical.pdf"),
  combined_two_panel_vertical,
  width = 18,
  height = 25,
  units = "in"
)
