BEGIN{
   print_on = 1;
}

($1=="/*__BEGIN_EXPERIMENTAL_SECTION__*/"){
   print_on = 0;
}

($1=="#__BEGIN_EXPERIMENTAL_SECTION__#"){
   print_on = 0;
}

($1=="/*___END_EXPERIMENTAL_SECTION___*/"){
   print_on = 1;
}

($1=="#___END_EXPERIMENTAL_SECTION___#"){
   print_on = 1;
}

($1=="/*UNCOMMENT"){
   getline;
   getline;
   while ($1 != "#endif"){
      print;
      getline;
   }
   getline;
}

($1=="/*" && $3=="Copyright"){
  getline;
  printf("/* (c) Copyright 2000-2003 Ted Ralphs. All Rights Reserved.                  */\n");
}

($1=="#" && $3=="Copyright"){
  getline;
  printf("# (c) Copyright 2000-2003 Ted Ralphs. All Rights Reserved.                   #\n");
}
($1!="/*___END_EXPERIMENTAL_SECTION___*/" && $1!="/*__BEGIN_EXPERIMENTAL_SECTION__*/" && $1!="/*UNCOMMENT*/" && $1!="#___END_EXPERIMENTAL_SECTION___#" &&
$1!="#___END_EXPERIMENTAL_SECTION___#"){
   if (print_on){
      print;
   }
}

END{
}







