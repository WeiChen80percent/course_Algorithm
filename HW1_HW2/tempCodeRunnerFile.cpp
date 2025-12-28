fprintf(gnuplotpipe, "set terminal pngcairo\n");
    fprintf(gnuplotpipe, "set title '%s'\n",output_file_address.c_str());
    fprintf(gnuplotpipe, "set output '%s.png'\n",output_file_address.c_str());
    fprintf(gnuplotpipe, "set xrange [0:100]\n");
    fprintf(gnuplotpipe, "set yrange [0:100]\n");
    fprintf(gnuplotpipe, "unset key\n");
    fprintf(gnuplotpipe, "plot '-' lt 1 lc 1 w lp\n");
    for(int i=0;i<min_nodeset.size();i++){
        fprintf(gnuplotpipe, "%d %d\n", min_nodeset[i].x_coord,min_nodeset[i].y_coord);
    }
    fprintf(gnuplotpipe, "e\n");