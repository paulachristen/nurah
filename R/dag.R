#' Define a Directed Acyclic Graph (DAG) for mortality prediction
#'
#' @param nodes A vector of node names (health domains).
#' @param edges A data frame specifying directed edges between nodes with columns from and to.
#' @param parameters A data frame specifying parameters for each edge including:
#'   - from: node from which effect originates
#'   - to: node affected
#'   - effect_size: numeric magnitude of effect
#'   - lag: numeric time lag in days (default 0)
#'
#' @return A list containing:
#'   - dag: A dagitty object representing the DAG.
#'   - parameters: parameters dataframe.
#'
#' @export
#'
#' @examples
#' nodes <- c("Nutrition", "Food insecurity", "Mortality")
#' edges <- data.frame(from = c("Food insecurity", "Nutrition"),
#'                     to = c("Nutrition", "Mortality"))
#' parameters <- data.frame(from = c("Food insecurity", "Nutrition"),
#'                          to = c("Nutrition", "Mortality"),
#'                          effect_size = c(-0.5, 1.2),
#'                          lag = c(30, 0))
#' dag <- define_dag(nodes, edges, parameters)
define_dag <- function(nodes, edges, parameters = NULL) {

  # create daggity strings
  edges_str <- apply(edges, 1, function(x) paste0(x["from"], " -> ", x["to"]))
  dagitty_str <- paste("dag {", paste(edges_str, collapse="; "), "}")

  # create daggity object and add parameters
  dag <- dagitty::dagitty(dagitty_str)
  dag_obj <- list(
    dag = dag,
    parameters = parameters
  )

  class(dag_obj) <- "nurah_dag"

  return(dag_obj)
}


#' Corrected Checchi et al. (2017) DAG Definition
#'
#' Defines the DAG from Checchi et al. (2017), accurately including all edges.
#'
#' @param parameters Optional dataframe with parameters for DAG edges, with columns:
#'   - from: origin node
#'   - to: destination node
#'   - effect_size: numeric magnitude of effect
#'   - lag: numeric time lag in days (default 0)
#'
#' @return A DAG based on Checchi et al. (2017) with class "nurah_dag".
#' @export
#'
#' @examples
#' params <- data.frame(
#'   from = c("Exposure to armed attacks or mechanical force of nature", "Food insecurity",
#'            "Nutritional status", "Burden of endemic infectious diseases"),
#'   to = c("Burden and typology of injuries", "Nutritional status",
#'          "Burden of endemic infectious diseases", "Population mortality"),
#'   effect_size = c(1.0, 0.7, 1.5, 2.0),
#'   lag = c(0, 30, 0, 0)
#' )
#' dag <- checchi_2017_dag(parameters = params)
checchi_2017_dag <- function(parameters = NULL) {

  nodes <- c(
    "\"Exposure to armed attacks or mechanical force of nature\"",
    "\"Forced displacement\"",
    "\"Interruption of chronic treatment\"",
    "\"Food insecurity\"",
    "\"Feeding and care practices\"",
    "\"Nutritional status\"",
    "\"Sexual and gender-based violence\"",
    "\"Mental health and psychosocial functioning\"",
    "\"Addiction\"",
    "\"Reproductive and neonatal health\"",
    "\"Burden of NCDs\"",
    "\"Burden of endemic infectious diseases\"",
    "\"Epidemic occurrence and severity\"",
    "\"Burden and typology of injuries\"",
    "\"Population mortality\"",
    "\"Humanitarian public health services\""
  )

  edges <- data.frame(
    from = c(
      "\"Exposure to armed attacks or mechanical force of nature\"", "\"Exposure to armed attacks or mechanical force of nature\"",
      "\"Exposure to armed attacks or mechanical force of nature\"", "\"Exposure to armed attacks or mechanical force of nature\"",
      "\"Exposure to armed attacks or mechanical force of nature\"", "\"Exposure to armed attacks or mechanical force of nature\"",
      "\"Forced displacement\"", "\"Forced displacement\"", "\"Forced displacement\"", "\"Forced displacement\"",
      "\"Forced displacement\"", "\"Food insecurity\"", "\"Food insecurity\"",
      "\"Feeding and care practices\"", "\"Interruption of chronic treatment\"", "\"Interruption of chronic treatment\"",
      "\"Sexual and gender-based violence\"", "\"Mental health and psychosocial functioning\"",
      "\"Nutritional status\"", "\"Nutritional status\"", "\"Nutritional status\"", "\"Nutritional status\"", "\"Nutritional status\"",
      "\"Addiction\"", "\"Burden of endemic infectious diseases\"", "\"Burden of NCDs\"",
      "\"Epidemic occurrence and severity\"", "\"Burden and typology of injuries\"", "\"Mental health and psychosocial functioning\"",
      "\"Humanitarian public health services\""
    ),
    to = c(
      "\"Forced displacement\"", "\"Interruption of chronic treatment\"",
      "\"Sexual and gender-based violence\"", "\"Burden and typology of injuries\"",
      "\"Mental health and psychosocial functioning\"", "\"Humanitarian public health services\"",
      "\"Interruption of chronic treatment\"", "\"Sexual and gender-based violence\"",
      "\"Food insecurity\"", "\"Feeding and care practices\"",
      "\"Burden of endemic infectious diseases\"", "\"Feeding and care practices\"", "\"Nutritional status\"",
      "\"Nutritional status\"", "\"Burden of NCDs\"", "\"Burden of endemic infectious diseases\"",
      "\"Mental health and psychosocial functioning\"", "\"Addiction\"",
      "\"Reproductive and neonatal health\"", "\"Burden of NCDs\"",
      "\"Burden of endemic infectious diseases\"", "\"Epidemic occurrence and severity\"", "\"Population mortality\"",
      "\"Burden of NCDs\"", "\"Population mortality\"", "\"Population mortality\"",
      "\"Population mortality\"", "\"Population mortality\"", "\"Population mortality\"",
      "\"Population mortality\""
    )
  )


  # create positional layout
  layout_matrix <- matrix(c(
    1, 8,   # Exposure
    1, 6,   # Forced displacement
    4, 6,   # Interruption
    2, 4,   # Food insecurity
    2, 5,   # Feeding
    4, 4,   # Nutritional status
    5, 7,   # Sexual violence
    7, 7,   # Mental health
    5, 5,   # Addiction
    6, 4,   # Neonatal health
    6, 3,   # NCDs
    6, 2,   # Infectious diseases
    6, 1,   # Epidemic severity
    6, 0,   # Injuries
    9, 2,   # Mortality
    5, 8    # Humanitarian services
  ), ncol=2, byrow=TRUE, dimnames = list(NULL, c("x", "y")))

  layout <- data.frame("name" = gsub("\"", "", nodes), layout_matrix)

  # define the dag and layout attr
  dag <- define_dag(nodes, edges, parameters)
  attr(dag, "layout") <- layout
  return(dag)
}


#' Generate Dummy Parameters for Checchi et al. (2017) DAG
#'
#' Returns a data frame of plausible effect sizes and lags for use with the Checchi et al. (2017) mortality DAG.
#' The parameters represent estimated contributions to crude mortality rate (CMR) in deaths per 10,000 people per day.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{from}{Name of the predictor node}
#'   \item{to}{Name of the outcome node}
#'   \item{effect_size}{Change in CMR (deaths per 10,000/day) per 1-unit increase in the predictor}
#'   \item{lag}{Lag in days before the effect manifests}
#' }
#'
#' @examples
#' params <- dummy_checchi_2017_parameters()
#' head(params)
#'
#' @export
dummy_checchi_2017_parameters <- function() {
  data.frame(
    from = c(
      "Exposure to armed attacks or mechanical force of nature",
      "Exposure to armed attacks or mechanical force of nature",
      "Exposure to armed attacks or mechanical force of nature",
      "Exposure to armed attacks or mechanical force of nature",
      "Exposure to armed attacks or mechanical force of nature",
      "Exposure to armed attacks or mechanical force of nature",
      "Forced displacement",
      "Forced displacement",
      "Forced displacement",
      "Forced displacement",
      "Forced displacement",
      "Food insecurity",
      "Food insecurity",
      "Feeding and care practices",
      "Interruption of chronic treatment",
      "Interruption of chronic treatment",
      "Sexual and gender-based violence",
      "Mental health and psychosocial functioning",
      "Nutritional status",
      "Nutritional status",
      "Nutritional status",
      "Nutritional status",
      "Nutritional status",
      "Addiction",
      "Burden of endemic infectious diseases",
      "Burden of NCDs",
      "Epidemic occurrence and severity",
      "Burden and typology of injuries",
      "Mental health and psychosocial functioning",
      "Humanitarian public health services"
    ),
    to = c(
      "Forced displacement",
      "Interruption of chronic treatment",
      "Sexual and gender-based violence",
      "Burden and typology of injuries",
      "Mental health and psychosocial functioning",
      "Humanitarian public health services",
      "Interruption of chronic treatment",
      "Sexual and gender-based violence",
      "Food insecurity",
      "Feeding and care practices",
      "Burden of endemic infectious diseases",
      "Feeding and care practices",
      "Nutritional status",
      "Nutritional status",
      "Burden of NCDs",
      "Burden of endemic infectious diseases",
      "Mental health and psychosocial functioning",
      "Addiction",
      "Reproductive and neonatal health",
      "Burden of NCDs",
      "Burden of endemic infectious diseases",
      "Epidemic occurrence and severity",
      "Population mortality",
      "Burden of NCDs",
      "Population mortality",
      "Population mortality",
      "Population mortality",
      "Population mortality",
      "Population mortality",
      "Population mortality"
    ),
    effect_size = c(
      0.2, 0.3, 0.4, 1.2, 0.2, -0.5,
      0.4, 0.4, 0.5, 0.4, 0.5,
      0.4, 0.6, 0.3,
      0.5, 0.6,
      0.4, 0.3,
      0.8, 0.7, 1.8, 0.6, 1.8,
      0.4, 2.0, 0.7, 2.5, 1.2, 0.3, -1.5
    ),
    lag = c(
      0, 14, 14, 0, 30, 0,
      14, 14, 30, 14, 45,
      0, 30, 0,
      30, 60,
      14, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0
    ),
    stringsAsFactors = FALSE
  )
}



#' Visualise a Directed Acyclic Graph (DAG) with Rounded‐Corner Boxes
#'
#' Plots a dag, using a predefined layout if available.
#'
#' @param dag A nurah_dag from `define_dag()` or `checchi_2017_dag()`.
#' @param use_dag_layout_attr Logical; if TRUE, use attr(dag, "layout") if present.
#' @param node_color Colour for node fill. Default \"skyblue\".
#' @param label_size Numeric cex for labels. Default 5
#'
#' @return Invisibly, the matrix of (x,y) coordinates used.
#' @export
visualise_dag <- function(dag,
                          use_dag_layout_attr = TRUE,
                          node_color = "skyblue",
                          label_size = 5) {

  if (!inherits(dag, "nurah_dag"))
    stop("Input must be a nurah_dag")

  # use default layout if present
  if (use_dag_layout_attr && !is.null(attr(dag, "layout"))) {
    tidy_dag <- tidy_dagitty_coords(dag$dag, coords = attr(dag, "layout"))
  } else {
    tidy_dag <- ggdag::tidy_dagitty(dag$dag)
  }

  # create default plot - use ggplot base to avoid double-drawing edges
  ggplot2::ggplot(tidy_dag, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
    ggdag::geom_dag_edges(edge_color = "gray40") +
    ggdag::geom_dag_point(color = node_color, size = 8) +
    ggdag::geom_dag_label(ggplot2::aes(label = .data$name), size = label_size) +
    ggdag::theme_dag() +
    ggplot2::labs(title = "Directed Acyclic Graph (DAG)")
}

#' Tidy a dagitty object with optional coordinates
#'
#' Converts a dagitty DAG into a tidy data frame. Allows user-provided coordinates,
#' otherwise defaults to automatic layout.
#'
#' @param .dagitty A dagitty object representing a DAG.
#' @param coords Optional. A named list of coordinates for nodes. If not provided, coordinates are generated.
#' @param seed Optional numeric seed for reproducible layout.
#' @param layout A layout type; defaults to `"nicely"`. `"time_ordered"` is also supported.
#' @param ... Additional arguments passed to the layout generation functions.
#'
#' @return A tidy dagitty object suitable for use with ggdag functions.
tidy_dagitty_coords <- function(.dagitty, coords = NULL, seed = NULL, layout = "nicely", ...) {
  if (!is.null(seed))
    set.seed(seed)

  if (dagitty::graphType(.dagitty) != "dag") {
    stop("`.dagitty` must be of graph type `dag`")
  }

  dag_edges <- ggdag:::get_dagitty_edges(.dagitty)

  if (!is.null(coords)) {
    # Use user-provided coordinates
    dagitty::coordinates(.dagitty) <- coords %>% ggdag:::coords2list()
  } else if (layout == "time_ordered") {
    coords <- dag_edges %>%
      ggdag:::edges2df() %>%
      ggdag:::auto_time_order() %>%
      ggdag:::time_ordered_coords() %>%
      ggdag:::coords2list()

    dagitty::coordinates(.dagitty) <- coords
  } else {
    ggdag:::check_verboten_layout(layout)
  }

  coords_df <- dag_edges %>%
    dplyr::select(.data$name, .data$to) %>%
    ggdag:::generate_layout(layout = layout,
                            vertices = names(.dagitty),
                            coords = dagitty::coordinates(.dagitty),
                            ...)

  tidy_dag <- dag_edges %>% ggdag:::tidy_dag_edges_and_coords(coords_df)

  ggdag:::new_tidy_dagitty(tidy_dag, .dagitty)
}


#' Topologically Sort DAG Nodes
#'
#' Returns a topologically sorted vector of node names from a directed acyclic graph
#' (DAG), ensuring that each node appears after all of its ancestor (parent) nodes.
#' This is useful for determining an order in which to simulate or process nodes so
#' that dependencies are respected.
#'
#' @param dag A DAG object of class \code{nurah_dag} (as created by \code{define_dag()})
#' containing the DAG structure and parameters.
#'
#' @return A character vector of node names sorted in topological order (parents always
#' come before their children).
#'
#' @details The function extracts the list of parent relationships from the \code{dag}
#' and then performs a topological sort using Kahn's algorithm. If the graph is not a
#' valid DAG (for example, if it contains a cycle or an undefined parent reference),
#' the function will throw an error.
#'
#' All parent names must correspond to nodes defined in the DAG. If any parent node is
#' referenced that is not in the DAG's node set, an error is produced.
#'
#' @examples
#' # Define a simple DAG: A -> B -> C
#' nodes <- c("A", "B", "C")
#' edges <- data.frame(from = c("A", "B"), to = c("B", "C"))
#' params <- data.frame(from = c("A", "B"), to = c("B", "C"), effect_size = c(0.5, 1.0), lag = c(0, 0))
#' dag_obj <- define_dag(nodes, edges, params)
#' topological_sort_dag(dag_obj)
#' #> [1] "A" "B" "C"
#' @export
topological_sort_dag <- function(dag) {
  if (!inherits(dag, "nurah_dag")) {
    stop("topological_sort_dag: Input must be a 'nurah_dag' object.")
  }
  # Extract all nodes and parent relationships from the dag
  params_df <- dag$parameters
  # Determine the set of all nodes (include nodes that might have no incoming or outgoing edges)
  if (!is.null(params_df) && nrow(params_df) > 0) {
    all_nodes <- unique(c(as.character(params_df$from), as.character(params_df$to)))
  } else {
    # If no edges defined in parameters, try to get nodes from dagitty object
    all_nodes <- tryCatch(dagitty::names(dag$dag), error = function(e) NULL)
    if (is.null(all_nodes)) {
      stop("DAG contains no nodes or edges to sort.")
    }
  }
  # Build parent list for each node
  parent_list <- setNames(vector("list", length(all_nodes)), all_nodes)
  if (!is.null(params_df) && nrow(params_df) > 0) {
    for (i in seq_len(nrow(params_df))) {
      from_node <- as.character(params_df$from[i])
      to_node   <- as.character(params_df$to[i])
      parent_list[[to_node]] <- c(parent_list[[to_node]], from_node)
    }
  }
  # Replace NULLs with character(0) for nodes with no parents
  parent_list <- lapply(parent_list, function(parents) {
    if (is.null(parents)) character(0) else as.character(parents)
  })

  # Verify that all parents are defined as nodes in the DAG
  all_parents <- unique(unlist(parent_list))
  undefined_parents <- setdiff(all_parents, all_nodes)
  if (length(undefined_parents) > 0) {
    stop("Undefined parent nodes in DAG: ", paste(undefined_parents, collapse = ", "))
  }

  # Perform Kahn's algorithm for topological sorting
  indegree <- setNames(integer(length(all_nodes)), all_nodes)
  # Calculate indegree for each node
  for (node in all_nodes) {
    for (parent in parent_list[[node]]) {
      indegree[node] <- indegree[node] + 1
    }
  }
  # Initialize list S with all nodes having no parents (indegree 0)
  S <- sort(all_nodes[indegree == 0])
  sorted <- character(0)

  while (length(S) > 0) {
    # Remove the first node from S
    n <- S[1]
    S <- S[-1]
    # Append it to the sorted list
    sorted <- c(sorted, n)
    # Decrease indegree of each child of n by 1
    for (child in all_nodes) {
      if (n %in% parent_list[[child]]) {
        indegree[child] <- indegree[child] - 1
        if (indegree[child] == 0) {
          # If child has no remaining parents, add to S
          S <- sort(c(S, child))
        }
      }
    }
  }

  if (length(sorted) != length(all_nodes)) {
    stop("DAG contains a cycle or undefined dependencies; topological sort failed.")
  }
  return(sorted)
}
