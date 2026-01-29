library(orderly)

orderly_run("1_data_cleaning")

orderly_run("2_model_fit")

orderly_run("3_plot")

orderly_validate_archive(action = "orphan")

orderly_prune_orphans()
