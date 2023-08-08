#!/bin/bash

script_name=$1

case "$script_name" in
    delete_tracks)
		echo 0
		;;

    unzip_data)
        echo 1
        ;;

    save_individual_images)
        echo 2
        ;;
     
    find_centerlines)
        echo 3
        ;;

    auto_resolve_problems)
        echo 4
        ;;

    calculate_spectra)
        echo 5
        ;;

    calculate_embeddings)
        echo 6
        ;;

	calculate_behaviors)
        echo 7
        ;;

    stimulus_realign)
        echo 8
        ;;

	plot_image_directory)
        echo 9
        ;;

    zip_data)
        echo 10
        ;;

    max)
        echo 10
        ;;
        
    *)
        echo -1
        ;;
esac

