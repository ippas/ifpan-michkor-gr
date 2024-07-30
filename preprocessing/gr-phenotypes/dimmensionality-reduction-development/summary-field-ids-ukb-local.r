install.packages("progress")
install.packages("furrr")
install.packages("writexl")

library(tidyverse)
library(magrittr)
library(httr)
library(rvest)
library(stringr)
library(progress)
library(purrr)
library(furrr)
library(writexl)

options(width = 1000)


# Define the function with retry logic
get_name_category <- function(category_id, attempts = 5, delay = 5) {
  url <- paste0("https://biobank.ctsu.ox.ac.uk/showcase/label.cgi?id=", category_id)

  for (i in 1:attempts) {
    try(
      {
        # Wczytaj zawartość strony
        page <- read_html(url)

        # Wyciągnij całą ścieżkę kategorii z div#main
        path_text <- page %>%
          html_node("div#main") %>%
          html_text()

        # Podziel tekst na kategorie używając strzałki jako separatora
        categories_all <- path_text %>%
          str_split("\\u9205") %>%
          unlist() %>%
          str_split("\n") %>%
          unlist() %>%
          .[1] %>%
          trimws() %>%
          str_replace_all("Category ", "") %>%
          str_replace_all(as.character(category_id), "") %>%
          str_split("⏵") %>%
          .[[1]] %>%
          str_trim()

        # get last element from vector categories_all
        last_category <- categories_all[length(categories_all)]

        return(last_category)
      },
      silent = TRUE
    )

    # If an error occurs, wait for 'delay' seconds before the next attempt
    Sys.sleep(delay)
  }

  # If all attempts fail, return NULL
  return(NULL)
}

################################################################################
# read data
################################################################################

# read table ukb-phenotypes-processed-raw2.tsv file to dataframe
summary_fields_ares <- read_tsv(
  "results/google-drive/preparing-phenotypes/ukb-phenotypes-processed-raw2.tsv"
)
summary_fields_ares

summary_fields_all_ukb <- read_tsv(
  "results/google-drive/preparing-phenotypes/field.tsv"
)

fields_categories <- read_tsv(
  "results/google-drive/preparing-phenotypes/fields_categories.tsv"
)

data_coding_ukb_metadata <- read_tsv(
  "results/google-drive/preparing-phenotypes/dimmensionality-reduction-development/data-coding-ukb.tsv"
)

data_coding_ukb_metadata %>%
  group_by(data_coding) %>%
  nest() %>%
  mutate(num_questions = map_int(data, ~ nrow(.x))) %>%
  select(!data) -> data_coding_ukb_metadata


################################################################################
# preprocessing data
################################################################################

# get unique main categories
summary_fields_all_ukb %>%
  select(main_category) %>%
  unique() %>%
  mutate(
    name_main_category = map(
      main_category,
      get_name_category
    )
  ) %>%
  unnest(name_main_category) %>%
  as.data.frame() -> map_name_main_category

# join fields_categories with map_name_main_category by name_main_category and last_category
fields_categories %>%
  left_join(
    .,
    map_name_main_category,
    by = c("last_category" = "name_main_category")
  ) %>%
  unique() -> fields_categories

fields_categories %>%
  na.omit() %>%
  dim()

left_join(
  summary_fields_all_ukb,
  map_name_main_category,
  by = "main_category"
) %>%
  mutate(
    on_ares = ifelse(
      field_id %in% as.character(summary_fields_ares$field_id),
      "yes",
      "no"
    )
  ) -> preprocessed_fields

left_join(
  summary_fields_all_ukb,
  map_name_main_category,
  by = "main_category"
) %>%
  mutate(
    on_ares = ifelse(
      field_id %in% as.character(summary_fields_ares$field_id),
      "yes",
      "no"
    )
  ) %>%
  .$main_category %>%
  unique() %>%
  length()

###############################################################################
# clean field IDs
###############################################################################
preprocessed_fields %>%
  select(on_ares) %>%
  table()

preprocessed_fields %>% dim()

# removing "firts reported" association field id, with field_id >= 1300000
preprocessed_fields %>%
  filter(field_id < 130000) -> preprocessed_fields

# removing redundant ICD9/10 categories
preprocessed_fields %>%
  # removing redundant ICD9 categories
  filter(!(main_category %in% c(
    40008, 40012, 40019, 40021, 40005, 4001, 40004, 40013
  ))) %>%
  # removing redundant ICD10 categories
  filter(!(main_category %in% c(
    40002, 41202, 41204, 41201, 40006, 40001, 41142, 41104, 41078
  ))) -> preprocessed_fields

preprocessed_fields %>% dim()

# removing technical categories
preprocessed_fields %>%
  filter(!(main_category %in% c(
    301, 302, 23165, 268, 265, 186, 187, 270, 185, 172, 171,
    100315, 100035, 100313, 199001, 100319,
    2000, 3001, 2006, 1020, 2005, 300, 170, 100069, 100070
  ))) -> preprocessed_fields

preprocessed_fields %>%
  select(on_ares) %>%
  table()

# removing redundant filed ids main category 2005 (Summary Operations), without field id equal 41272
preprocessed_fields %>%
  filter(!(main_category == 2005 & field_id != 41272)) -> preprocessed_fields

# removing 2006 and 3001 main categories (technical categories)
preprocessed_fields %>%
  filter(!(main_category %in% c(2006, 3001))) -> preprocessed_fields

# removing categories associated with local environment
preprocessed_fields %>%
  filter(!(main_category %in% c(150, 114, 115, 151, 603))) -> preprocessed_fields

# removing  private fields
preprocessed_fields %>%
  filter(private != 1) -> preprocessed_fields

# removing field ids with num_participants < 50000
preprocessed_fields %>%
  filter(num_participants >= 50000) %>%
  select(on_ares) %>%
  table()

preprocessed_fields %>%
  filter(num_participants >= 50000) -> preprocessed_fields

preprocessed_fields %>% dim()

preprocessed_fields %>%
  filter(num_participants > 490000) %>%
  dim()

preprocessed_fields %>%
  filter(num_participants <= 490000) -> preprocessed_fields

preprocessed_fields %>%
  select(on_ares) %>%
  table()

# prepare function to save table to tsv and xlsx file
save_tsv_xlsx <- function(df, file_name) {
  df %>%
    write_tsv(
      paste0(file_name, ".tsv")
    )

  df %>%
    write_xlsx(
      paste0(file_name, ".xlsx")
    )
}

preprocessed_fields %>%
  as.data.frame() %>%
  select(!c(
    "availability", "stability", "private", "value_type",
    "base_type", "item_type", "strata", "arrayed", "sexed", "notes",
    "cost_do", "cost_on", "cost_sc", "array_min", "array_max",
    "instanced", "instance_id", "instance_min", "instance_max"
  )) %>%
  select(c(
    "field_id", "title", "main_category", "name_main_category",
    "encoding_id", "units", "num_participants", "item_count",
    "showcase_order", "on_ares", "debut", "version"
  )) %>%
  save_tsv_xlsx(
    ., "results/google-drive/preparing-phenotypes/dimmensionality-reduction-development/summary-fields-ukb"
  )


preprocessed_fields %>%
  as.data.frame() %>%
  select(!c(
    "availability", "stability", "private", "value_type",
    "base_type", "item_type", "strata", "arrayed", "sexed", "notes",
    "cost_do", "cost_on", "cost_sc", "array_min", "array_max",
    "instanced", "instance_id", "instance_min", "instance_max"
  )) %>%
  select(c(
    "field_id", "title", "main_category", "name_main_category",
    "encoding_id", "units", "num_participants", "item_count",
    "showcase_order", "on_ares", "debut", "version"
  )) %>%
  left_join(
    ., data_coding_ukb_metadata,
    by = c("encoding_id" = "data_coding")
  ) %>%
  # in column num_question replace NA with 1
  mutate(num_questions = ifelse(is.na(num_questions), 1, num_questions)) %>%
  select(!c(debut, version, on_ares, units, title, showcase_order)) %>%
  group_by(main_category, name_main_category) %>%
  nest() %>%
  mutate(sum_field = map_int(data, ~ nrow(.x))) %>%
  mutate(sum_questions = map_int(data, ~ sum(.x$num_questions))) %>%
  mutate(min_parcitipants = map_dbl(data, ~ min(.x$num_participants))) %>%
  mutate(max_parcitipants = map_dbl(data, ~ max(.x$num_participants))) %>%
  mutate(range_participants = max_parcitipants - min_parcitipants) %>%
  select(!data) %>%
  save_tsv_xlsx(
    ., "results/google-drive/preparing-phenotypes/dimmensionality-reduction-development/summary-category-ukb-metadata"
  )



summary_fields_all_ukb$encoding_id %>%
  unique() %>%
  length()
# select(c(field_id, main_category, name_main_category, num_participants, on_ares)) %>%
# write_tsv(
#   "results/google-drive/preparing-phenotypes/summary-fields-ukb.tsv"
# )


################################################################
# testing code
################################################################

preprocessed_fields %>%
  select(c(field_id, title, main_category, name_main_category, num_participants, on_ares)) %>%
  group_by(main_category) %>%
  nest() %>%
  mutate(max_parcitipants = map_dbl(data, ~ max(.x$num_participants))) %>%
  mutate(min_parcitipants = map_dbl(data, ~ min(.x$num_participants))) %>%
  mutate(range = max_parcitipants - min_parcitipants) %>%
  filter(range < 100) %>%
  unnest(data) %>%
  mutate(max = ifelse(num_participants == max_parcitipants, "yes", "no")) %>%
  filter(max == "yes") %>%
  nest() %>%
  mutate(data = map(data, ~ slice(.x, 1))) %>%
  unnest(data) %>%
  .$field_id -> fields_selected_by_range

preprocessed_fields %>%
  select(c(field_id, title, main_category, name_main_category, num_participants, on_ares)) %>%
  group_by(main_category) %>%
  nest() %>%
  mutate(max_parcitipants = map_dbl(data, ~ max(.x$num_participants))) %>%
  mutate(min_parcitipants = map_dbl(data, ~ min(.x$num_participants))) %>%
  mutate(range = max_parcitipants - min_parcitipants) %>%
  filter(range < 100) %>%
  unnest(data) %>%
  filter(!(field_id %in% fields_selected_by_range)) %>%
  .$field_id -> fileds_to_remove_by_range

preprocessed_fields %>%
  filter(!(field_id %in% fileds_to_remove_by_range)) %>%
  select(on_ares) %>%
  table()
