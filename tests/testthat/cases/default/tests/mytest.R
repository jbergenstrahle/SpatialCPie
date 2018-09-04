app <- ShinyDriver$new("../")
app$snapshotInit("mytest")

# Input 'tree_selected' was set, but doesn't have an input binding.
app$snapshot()
