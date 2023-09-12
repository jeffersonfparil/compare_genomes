workflow {
    list_files = file(projectDir + "/../work/*/*/.command.err", hidden: true)
    // println(list_files)
    // println(list_files.size())
    // for (f : list_files) {
    //     print f.text
    // }
    // Channel
    //     .from(list_files)
    //     .toList()
    //     .sort()
    //     .last()
    //     .view()
    print list_files[list_files.size()-1].text
}