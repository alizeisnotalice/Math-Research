ClearAll[h,u,w];
Print["Wolfram: ",$Version];
inc=Min[h,u+w]-Min[h,u];
Print["domain h,u,w>=0"];
Print["increment bounds violation: ",Reduce[h>=0&&u>=0&&w>=0&&(inc<0||inc>w),{h,u,w},Reals]];
Print["unsaturated equality violation: ",Reduce[h>=0&&u>=0&&w>=0&&h>=u+w&&inc!=w,{h,u,w},Reals]];
