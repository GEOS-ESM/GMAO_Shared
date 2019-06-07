! test_utils.f90
! Test utilities module

program test_utils
    use ncdw_strarrutils
    
    character(len=100) :: str1, str2
    character(len=:), allocatable :: strpart
    
    str1 = "Hello, world! Today's a great day!"
    str2 = "aabbababababbbbaaaaa!"
    
    print *, "String #1:"
    print *, str1
    
    print *, "String part before spaces:"
    call string_before_delimiter(trim(str1), " ", strpart)
    write (*, "(A, A, A)") "'", strpart, "'"
    
    print *, "String #2:"
    print *, str2
    print *, "Number of occurances of 'a':"
    print *, string_count_substr(str2, "a")
    
    print *, "Max len when split using 'a':"
    print *, string_get_max_split(str2, "a")
    
    print *, "Split using 'a':"
    call string_array_dump(string_split_index(str2, "a"))
    
    print *, "Number of occurances of 'aa':"
    print *, string_count_substr(str2, "aa")
    
    print *, "Max len when split using 'aa':"
    print *, string_get_max_split(str2, "aa")
    
    print *, "Split using 'aa':"
    call string_array_dump(string_split_index(str2, "aa"))
    
    print *, "Number of occurances of 'bb':"
    print *, string_count_substr(str2, "bb")
    
    print *, "Max len when split using 'bb':"
    print *, string_get_max_split(str2, "bb")
    
    print *, "Split using 'bb':"
    call string_array_dump(string_split_index(str2, "bb"))
    
    print *, "Number of occurances of 'ab':"
    print *, string_count_substr(str2, "ab")
    
    print *, "Max len when split using 'ab':"
    print *, string_get_max_split(str2, "ab")
    
    print *, "Split using 'ab':"
    call string_array_dump(string_split_index(str2, "ab"))
    
    print *, "Number of occurances of 'a!':"
    print *, string_count_substr(str2, "a!")
    
    print *, "Max len when split using 'a!':"
    print *, string_get_max_split(str2, "a!")
    
    print *, "Split using 'a!':"
    call string_array_dump(string_split_index(str2, "a!"))
end program test_utils
