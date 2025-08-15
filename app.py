import os
import matplotlib.pyplot as plt
import gradio as gr
import math
import ezdxf
import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import ttk
from io import BytesIO


def draw_quadrilateral_close_ended(area_length, area_width, cell_length, cell_width, baffle_width):
    # Convert grid dimensions from mm to meters
    cell_length = cell_length + baffle_width
    cell_width = cell_width + baffle_width
    
    cell_length_m = cell_length / 1000  # mm to meters
    cell_width_m = cell_width / 1000    # mm to meters

    # Initialize output variables
    quadrilaterals_close_grid_details = ""
    quadrilaterals_close_image_path = ""
    quadrilaterals_close_pdf_path = ""

    # Validate that grid dimensions do not exceed area dimensions
    if cell_length_m > area_length or cell_width_m > area_width:
       return "Error: Grid dimensions must be less than or equal to area dimensions.", None, None

    # Validate that baffle width does not exceed grid dimensions
    if baffle_width >= cell_length or baffle_width >= cell_width:
        return "Baffle width must be less than both Cell length and width.", None, None

    # Calculate number of main and cross members
    num_main_members = int(area_width / cell_width_m) + 1  # Horizontal lines
    num_cross_members = int(area_length / cell_length_m) + 1  # Vertical lines

    # Calculate initial free lengths
    free_length_main_member = (area_length - (num_cross_members - 1) * cell_length_m) / 2
    free_length_cross_member = (area_width - (num_main_members - 1) * cell_width_m) / 2

    # Adjust number of cross members based on free length of main members
    free_length_old_main = free_length_main_member

    if 0 <= free_length_old_main < 0.25 * cell_length_m:
        num_cross_members -= 2  # Reduce by two
        print("Adjusted number of cross members: Reduced by one due to free length condition.")

    elif 0.25 * cell_length_m <= free_length_old_main < 0.75 * cell_length_m:
        num_cross_members -= 1  # Reduce by two
        print("Adjusted number of cross members: Reduced by one due to free length condition.")

    elif 0.75 * cell_length_m <= free_length_old_main < 1.00 * cell_length_m:
        pass  # Do nothing

    # Recalculate the free lengths with the adjusted number of cross members
    free_length_main_member = (area_length - (num_cross_members - 1) * cell_length_m) / 2

    # Adjust number of main members based on free length of cross members
    free_length_old_cross = free_length_cross_member

    if 0 <= free_length_old_cross < 0.25 * cell_width_m:
        num_main_members -= 2
        print("Adjusted number of main members: Reduced by one due to cross member condition.")

    elif 0.25 * cell_width_m <= free_length_old_cross < 0.75 * cell_width_m:
        num_main_members -= 1
        print("Adjusted number of main members: Reduced by one due to cross member condition.")

    elif 0.75 * cell_width_m <= free_length_old_cross < 1.00 * cell_width_m:
        pass
    
    # Recalculate the free lengths with the adjusted number of main members
    free_length_cross_member = (area_width - (num_main_members - 1) * cell_width_m) / 2
    #cross_member_df = pd.DataFrame(columns=["Product", "Quantity (NoS)"])

    i = 0
    standard_length = 0
    while True:
        next_length = cell_length_m * (i + 1) * 1000
        if next_length > 3000:
            break
        
        i += 1
        standard_length = cell_length_m * i * 1000
    
    n = int(area_length/3)
    main_member_df = pd.DataFrame(columns=[
        "Main Member line", "Total Length (in mm)"
        ] + [
            item for k in range(n+1) for item in [
                f"Member Name {k+1}",
                f"Length {k+1} (in mm)",
                f"First Punch {k+1} (in mm)"
            ]
        ])
    i = 1
    total_main_member = 0
    for i in range(num_main_members):
        new_row = pd.DataFrame(
            [[f"M{i + 1}", f"{area_length*1000-baffle_width:.2f}"]],
            columns=["Main Member line", "Total Length (in mm)"]
        )
        main_member_df = pd.concat([main_member_df, new_row], ignore_index=True)
        total_length = area_length * 1000 - baffle_width
        solution_found = False
        tolerance = 0.01  # Acceptable rounding error  

        if total_length < 3000:                    
            main_member_df.at[main_member_df.index[-1], f"Member Name {1}"] = f"M{i + 1}_1"
            main_member_df.at[main_member_df.index[-1], f"Length {1} (in mm)"] = f"{total_length:.2f}"
            main_member_df.at[main_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(0-baffle_width/2):.2f}"
        else:
            A = cell_length_m * 500   
            B = free_length_main_member*1000 + cell_length_m * 500 - baffle_width/2

            y_max = int((3000 - B) // (cell_length_m * 1000))
            
            for y in reversed(range(y_max + 1)):
                first_part_length = B + y * cell_length_m * 1000
                print("First part Length", first_part_length)
                z_max = int((3000 - B) // (cell_length_m * 1000))
                for z in reversed(range(z_max + 1)):
                    last_part_length = B + z * cell_length_m * 1000
                    print ("Last part length = ", last_part_length)
                    remaining = total_length - first_part_length - last_part_length
                    n = remaining / standard_length
                    
                    if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                        n = round(n)
                        total_check = first_part_length + n * standard_length + last_part_length
                        if abs(total_length - total_check) <= tolerance:
                            main_member_df.at[main_member_df.index[-1], f"Member Name {1}"] = f"M{i + 1}_1"
                            main_member_df.at[main_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                            main_member_df.at[main_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{((cell_length-baffle_width)/2):.2f}"
                            
                            # Add last part (comes after n standard parts + first)
                            main_member_df.at[main_member_df.index[-1], f"Member Name {n+2}"] = f"M{i + 1}_{n+2}"
                            main_member_df.at[main_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                            main_member_df.at[main_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{((cell_length-baffle_width)/2):.2f}"      
                            for q in range(n):
                                main_member_df.at[main_member_df.index[-1], f"Member Name {q+2}"] = f"M{i + 1}_{q+2}"
                                main_member_df.at[main_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                main_member_df.at[main_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{((cell_length-baffle_width)/ 2):.2f}"                                      
                            total_main_member = total_main_member + n + 2
                            solution_found = True                                    
                            break

                if solution_found:
                    break

            if not solution_found:
                print("No valid combination found.")
                
	
	
    standard_length = 0
    i = 0
    while True:
        next_length = cell_width_m * (i + 1) * 1000
        if next_length > 3000:
            break
        
        i += 1
        standard_length = cell_width_m * i * 1000
    #print (" standard length= ", standard_length)
    
    n = int(area_width/3)
    cross_member_df = pd.DataFrame(columns=[
        "Cross Member line", "Total Length (in mm)"
        ] + [
            item for k in range(n+1) for item in [
                f"Member Name {k+1}",
                f"Length {k+1} (in mm)",
                f"First Punch {k+1} (in mm)"
            ]
        ])
    i = 1
    total_cross_member = 0
    for i in range(num_cross_members):
        new_row = pd.DataFrame(
            [[f"C{i + 1}", f"{area_width*1000 - baffle_width:.2f}"]],
            columns=["Cross Member line", "Total Length (in mm)"]
        )
        cross_member_df = pd.concat([cross_member_df, new_row], ignore_index=True)
        total_length = area_width * 1000 - baffle_width
        #print("total length = ",  total_length)
        solution_found = False
        tolerance = 0.01  # Acceptable rounding error  

        if total_length < 3000:                    
            cross_member_df.at[cross_member_df.index[-1], f"Member Name {1}"] = f"C{i + 1}_1"
            cross_member_df.at[cross_member_df.index[-1], f"Length {1} (in mm)"] = f"{total_length:.2f}"
            cross_member_df.at[cross_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(0-baffle_width):.2f}"
        else:
            A = cell_width_m * 500
            B = free_length_cross_member*1000 + cell_width_m * 500 - baffle_width/2

            y_max = int((3000 - B) // (cell_width_m * 1000))
            
            for y in reversed(range(y_max + 1)):
                first_part_length = B + y * cell_width_m * 1000 
                #print("first_part_length = ",  first_part_length)
                z_max = int((3000 - B) // (cell_width_m * 1000))
                for z in reversed(range(z_max + 1)):
                    last_part_length = B + z * cell_width_m * 1000
                    #print("last_part_length = ",  last_part_length) 
                    remaining = total_length - first_part_length - last_part_length
                    n = remaining / standard_length
                    
                    if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                        n = round(n)
                        total_check = first_part_length + n * standard_length + last_part_length
                        if abs(total_length - total_check) <= tolerance:
                            cross_member_df.at[cross_member_df.index[-1], f"Member Name {1}"] = f"C{i + 1}_1"
                            cross_member_df.at[cross_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                            cross_member_df.at[cross_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{((cell_width - baffle_width)/2):.2f}"
                            
                            # Add last part (comes after n standard parts + first)
                            cross_member_df.at[cross_member_df.index[-1], f"Member Name {n+2}"] = f"C{i + 1}_{n+2}"
                            cross_member_df.at[cross_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                            cross_member_df.at[cross_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{((cell_width - baffle_width)/2):.2f}"      
                            for q in range(n):
                                cross_member_df.at[cross_member_df.index[-1], f"Member Name {q+2}"] = f"C{i + 1}_{q+2}"
                                cross_member_df.at[cross_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                cross_member_df.at[cross_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{((cell_width - baffle_width)/2):.2f}"                                      
                            total_cross_member = total_cross_member + n + 2
                            solution_found = True                                    
                            break

                if solution_found:
                    break

            if not solution_found:
                print("No valid combination found.")

    i = 0
    num_outer_members = 2
    standard_length_m = 0
    while True:
        next_length = cell_length_m * (i + 1) * 1000
        if next_length > 3000:
            break
        
        i += 1
        standard_length_m = cell_length_m * i * 1000
    
    n = int(max(area_length, area_width)/3)
    outer_member_df = pd.DataFrame(columns=[
        "Outer Member line", "Total Length (in mm)"
        ] + [
            item for k in range(n+1) for item in [
                f"Member Name {k+1}",
                f"Length {k+1} (in mm)",
                f"First Punch {k+1} (in mm)"
            ]
        ])
    i = 1
    total_outer_member = 0
    for i in range(num_outer_members):
        new_row = pd.DataFrame(
            [[f"O{i + 1}", f"{area_length*1000:.2f}"]],
            columns=["Outer Member line", "Total Length (in mm)"]
        )
        outer_member_df = pd.concat([outer_member_df, new_row], ignore_index=True)
        total_length = area_length * 1000
        solution_found = False
        tolerance = 0.01  # Acceptable rounding error  

        if total_length < 3000:                    
            outer_member_df.at[outer_member_df.index[-1], f"Member Name {1}"] = f"O{i + 1}_1"
            outer_member_df.at[outer_member_df.index[-1], f"Length {1} (in mm)"] = f"{total_length:.2f}"
            outer_member_df.at[outer_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(free_length_main_member*1000-baffle_width):.2f}"
        else:
            A = cell_length_m * 500
            B = free_length_main_member*1000 + A

            y_max = int((3000 - A - B) // (cell_length_m * 1000))
            
            for y in reversed(range(y_max + 1)):
                first_part_length = A + B + y * cell_length_m * 1000 
                #print("First part Length", first_part_length)
                z_max = int((3000 - A - B) // (cell_length_m * 1000))
                for z in reversed(range(z_max + 1)):
                    last_part_length = A + B + z * cell_length_m * 1000
                    #print ("Last part length = ", last_part_length)  
                    remaining = total_length - first_part_length - last_part_length
                    n = remaining / standard_length_m
                    
                    if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                        n = round(n)
                        total_check = first_part_length + n * standard_length_m + last_part_length
                        if abs(total_length - total_check) <= tolerance:
                            outer_member_df.at[outer_member_df.index[-1], f"Member Name {1}"] = f"O{i + 1}_1"
                            outer_member_df.at[outer_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                            outer_member_df.at[outer_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{((cell_length - baffle_width)/2):.2f}"
                            
                            # Add last part (comes after n standard parts + first)
                            outer_member_df.at[outer_member_df.index[-1], f"Member Name {n+2}"] = f"O{i + 1}_{n+2}"
                            outer_member_df.at[outer_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                            outer_member_df.at[outer_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{((cell_length - baffle_width)/2):.2f}"      
                            for q in range(n):
                                outer_member_df.at[outer_member_df.index[-1], f"Member Name {q+2}"] = f"O{i + 1}_{q+2}"
                                outer_member_df.at[outer_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length_m:.2f}"
                                outer_member_df.at[outer_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{((cell_length - baffle_width)/ 2):.2f}"                                      
                            total_outer_member = total_outer_member + n + 2
                            solution_found = True                                    
                            break

                if solution_found:
                    break

            if not solution_found:
                print("No valid combination found.")
                
	
	
    standard_length = 0
    i = 0
    while True:
        next_length = cell_width_m * (i + 1) * 1000
        if next_length > 3000:
            break
        
        i += 1
        standard_length_c = cell_width_m * i * 1000
    #print (" standard length= ", standard_length_c)

    i = 1
    for i in range(num_outer_members):
        new_row = pd.DataFrame(
            [[f"O{i + 3}", f"{area_width*1000:.2f}"]],
            columns=["Outer Member line", "Total Length (in mm)"]
        )
        outer_member_df = pd.concat([outer_member_df, new_row], ignore_index=True)
        total_length = area_width * 1000
        solution_found = False
        tolerance = 0.01  # Acceptable rounding error  

        if total_length < 3000:                    
            outer_member_df.at[outer_member_df.index[-1], f"Member Name {1}"] = f"O{i + 3}_1"
            outer_member_df.at[outer_member_df.index[-1], f"Length {1} (in mm)"] = f"{total_length:.2f}"
            outer_member_df.at[outer_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(free_length_cross_member*1000-baffle_width):.2f}"
        else:
            A = cell_width_m * 500
            B = free_length_cross_member*1000 + A

            y_max = int((3000 - A - B) // (cell_width_m * 1000))
            
            for y in reversed(range(y_max + 1)):
                first_part_length = A + B + y * cell_width_m * 1000 
                z_max = int((3000 - A - B) // (cell_width_m * 1000))
                for z in reversed(range(z_max + 1)):
                    last_part_length = A + B + z * cell_width_m * 1000
                    #print ("Last part length = ", last_part_length)  
                    remaining = total_length - first_part_length - last_part_length
                    n = remaining / standard_length_c
                    
                    if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                        n = round(n)
                        total_check = first_part_length + n * standard_length_c + last_part_length
                        if abs(total_length - total_check) <= tolerance:
                            outer_member_df.at[outer_member_df.index[-1], f"Member Name {1}"] = f"O{i + 3}_1"
                            outer_member_df.at[outer_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                            outer_member_df.at[outer_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{((cell_width - baffle_width)/2):.2f}"
                            
                            # Add last part (comes after n standard parts + first)
                            outer_member_df.at[outer_member_df.index[-1], f"Member Name {n+2}"] = f"O{i + 3}_{n+2}"
                            outer_member_df.at[outer_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                            outer_member_df.at[outer_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{((cell_width - baffle_width)/2):.2f}"      
                            for q in range(n):
                                outer_member_df.at[outer_member_df.index[-1], f"Member Name {q+2}"] = f"O{i + 3}_{q+2}"
                                outer_member_df.at[outer_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length_c:.2f}"
                                outer_member_df.at[outer_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{((cell_width - baffle_width)/ 2):.2f}"                                      
                            total_outer_member = total_outer_member + n + 2
                            solution_found = True                                    
                            break

                if solution_found:
                    break

            if not solution_found:
                print("No valid combination found.")
    

 

    # Calculate total lengths after adjustments
    total_main_length = num_main_members * (area_length - baffle_width/1000)
    total_cross_length = num_cross_members * (area_width - baffle_width/1000)

    # Calculate running meter per square meter
    area = area_length * area_width
    running_meter_per_sqm = (total_main_length + total_cross_length + 2 * (area_length + area_width))  / area

      # Generate grid details text
    quadrilaterals_close_grid_details = f"""
    **Grid Details:**
    - Total Number of Main Member Lines: {num_main_members}
    - Total Number of Cross Member Lines: {num_cross_members}
    - Free Length of Main Member: {free_length_main_member * 1000:.2f} mm
    - Free Length of Cross Member: {free_length_cross_member * 1000:.2f} mm
    - Grid Size Opening: {cell_length - baffle_width} mm × {cell_width - baffle_width} mm
    - Total Length of Main Members: {total_main_length:.2f} meters
    - Total Length of Cross Members: {total_cross_length:.2f} meters
    - Total Length of Outer Members: {2 * (area_length + area_width):.2f} meters
    - Running Meter per Square Meter: {running_meter_per_sqm:.2f} running meters per square meter
    """

#def plot_grid(num_main_members, num_cross_members, area_width, area_length, free_length_cross_member, free_length_main_member, cell_width_m, cell_length_m):
    fig, ax = plt.subplots(figsize=(8.27, 11.69))
    doc = ezdxf.new()
    msp = doc.modelspace()

    i = 0

    for i in range(num_main_members):
        y = free_length_cross_member + i * cell_width_m
        ax.plot([0, area_length], [y, y], color='orange', label='Main Member' if i == 0 else "")
        msp.add_line((0, y*1000), (area_length*1000, y*1000    ))  # Save line to DXF
    j = 0 

    for j in range(num_cross_members):
        x = free_length_main_member + j * cell_length_m
        ax.plot([x, x], [0, area_width], color='green', label='Cross Member' if j == 0 else "")
        msp.add_line((x*1000, 0), (x*1000, area_width*1000    ))  # Save line to DXF

        # Plotting the additional lines
    ax.plot([0, area_length], [0, 0], color='blue', label='Outer Member')  # y = 0 line
    ax.plot([0, area_length], [area_width, area_width], color='blue', label='y=area_width' if num_main_members == 0 and num_cross_members == 0 else "")  # y = area_width line
    ax.plot([0, 0], [0, area_width], color='blue', label='x=0' if num_main_members == 0 and num_cross_members == 0 else "")  # x = 0 line
    ax.plot([area_length, area_length], [0, area_width], color='blue', label='x=area_length' if num_main_members == 0 and num_cross_members == 0 else "")  # x = area_length line

    # Create the outer members (DXF data)
    msp.add_line((0, 0), (area_length*1000, 0))  # x=0 line
    msp.add_line((0, area_width*1000), (area_length*1000, area_width*1000))  # x=area_length line
    msp.add_line((0, 0), (0, area_width*1000))  # y=0 line
    msp.add_line((area_length*1000, 0), (area_length*1000, area_width*1000))  # y=area_width line

    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(0, area_length + 1)
    ax.set_ylim(0, area_width + 1)
    ax.set_title(f"Grid Pattern of {cell_length} mm × {cell_width} mm ")
    ax.set_xlabel(f"Length of Area (m): {area_length} m")
    ax.set_ylabel(f"Width of Area (m): {area_width} m")

    ax.set_xticks([])  # Hide x-axis ticks
    ax.set_yticks([])  # Hide y-axis ticks
    # Hide the spines (the outer box)
    for spine in ax.spines.values():
        spine.set_visible(False)
    
    plt.grid(False)        # Disables grid lines
    legend_width = 0.5 * area_width
    legend_height = 0.5 * area_length
    # Transform data coordinates to figure coordinates
    transform = ax.transData.transform
    inv_transform = ax.transAxes.inverted().transform

    # Convert (area_length, area_width) from data coordinates to figure coordinates
    legend_x, legend_y = inv_transform(transform((area_length, area_width)))
    # Add legend with bottom-left corner at (area_length, area_width)
    ax.legend(loc='lower left', bbox_to_anchor=(legend_x, legend_y), bbox_transform=ax.transAxes, fontsize=8, frameon=True)
    plt.grid()
    doc.saveas("grid_plot.dxf")
    plt.show()

    # Define a safer directory path, e.g., within the current working directory
    output_dir = "./output"

    # Ensure the directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Set the full paths for saving
    quadrilaterals_close_pdf_path = os.path.join(output_dir, "grid_plot.pdf")
    quadrilaterals_close_image_path = os.path.join(output_dir, "grid_plot.png")
    quadrilaterals_close_dxf_path = os.path.join(output_dir, "grid_plot.dxf")
    

    # Save the figure as a PDF and PNG
    plt.savefig(quadrilaterals_close_pdf_path, format="pdf", bbox_inches="tight")  # Save as PDF
    plt.savefig(quadrilaterals_close_image_path, format="png", bbox_inches="tight")  # Save as PNG
    doc.saveas(quadrilaterals_close_dxf_path)  # Save DXF file
    plt.show()
    plt.close(fig)

    excel_output = BytesIO()
    with pd.ExcelWriter(excel_output, engine="openpyxl") as writer:
        main_member_df.to_excel(writer, index=False, sheet_name=f"Main Member- {total_main_member} nos")
        cross_member_df.to_excel(writer, index=False, sheet_name=f"Cross Members- {total_cross_member} nos")
        outer_member_df.to_excel(writer, index=False, sheet_name=f"Cross Members- {total_outer_member} nos")
    excel_output.seek(0)

    main_member_df = main_member_df.fillna('')
    cross_member_df = cross_member_df.fillna('')
    outer_member_df = outer_member_df.fillna('')

    # Generate HTML tables for display
    html_output = f"""
    <div style='font-family: Arial, sans-serif;'>    
        <h2>Main Members Details: {total_main_member} nos</h2>
        {main_member_df.to_html(index=False, classes='dataframe')}

        <h2>Cross Members Details: {total_cross_member} nos</h2>
        {cross_member_df.to_html(index=False, classes='dataframe')}

        <h2>Cross Members Details: {total_outer_member} nos</h2>
        {outer_member_df.to_html(index=False, classes='dataframe')}
        
        <style>
            .dataframe {{
                width: 100%;
                border-collapse: collapse;
                margin-bottom: 20px;
                color: #000000;
                background-color: #ffffff;
            }}
            .dataframe th {{
                background-color: #ffa8B6;
                color: #373435;
                padding: 10px;
                text-align: center;
            }}
            .dataframe td {{
                padding: 8px;
                border-bottom: 1px solid #ddd;
                text-align: center;
            }}
            .dataframe tr:nth-child(even) {{
                background-color: #f2f2f2;
            }}

            /* Dark theme overrides */
            @media (prefers-color-scheme: dark) {{
                .dataframe {{
                    color: #ffffff;
                    background-color: #121212;
                }}
                .dataframe th {{
                    background-color: #444;
                    color: #ffa8b6;
                }}
                .dataframe td {{
                    border-bottom: 1px solid #444;
                }}
                .dataframe tr:nth-child(even) {{
                    background-color: #1e1e1e;
                }}
            }}
        </style>
    </div>
    """

    # Save to disk for Gradio
    excel_file_path = os.path.join(output_dir, "triangular_grid_output.xlsx")
    with open(excel_file_path, "wb") as f:
        f.write(excel_output.getvalue())
    
    # ✅ Final return
    return quadrilaterals_close_grid_details, quadrilaterals_close_image_path, html_output, excel_file_path, quadrilaterals_close_pdf_path, quadrilaterals_close_dxf_path

def draw_quadrilaterals_open_ended(area_length, area_width, cell_length, cell_width, baffle_width):
    # Convert grid dimensions from mm to meters
    cell_length = cell_length + baffle_width
    cell_width = cell_width + baffle_width
    cell_length_m = cell_length / 1000  # mm to meters
    cell_width_m = cell_width / 1000    # mm to meters

    # Initialize output variables
    quadrilaterals_open_grid_details = ""
    quadrilaterals_open_image_path = ""
    quadrilaterals_open_pdf_path = ""

    # Validate that grid dimensions do not exceed area dimensions
    if cell_length_m > area_length or cell_width_m > area_width:
       return "Error: Grid dimensions must be less than or equal to area dimensions.", None, None

    # Validate that baffle width does not exceed grid dimensions
    if baffle_width >= cell_length or baffle_width >= cell_width:
        return "Baffle width must be less than both Cell length and width.", None, None

    # Calculate number of main and cross members
    num_main_members = int(area_width / cell_width_m) + 1  # Horizontal lines
    num_cross_members = int(area_length / cell_length_m) + 1  # Vertical lines

    # Calculate initial free lengths
    free_length_main_member = (area_length - (num_cross_members - 1) * cell_length_m) / 2
    free_length_cross_member = (area_width - (num_main_members - 1) * cell_width_m) / 2

    # Adjust number of cross members based on free length of main members
    free_length_old_main = free_length_main_member

    if 0 <= free_length_old_main < 0.10 * cell_length_m:
        num_cross_members -= 1  # Reduce by one
        print("Adjusted number of cross members: Reduced by one due to free length condition.")

    elif 0.10 * cell_length_m <= free_length_old_main < 0.20 * cell_length_m:
        if 0.60 * cell_length_m > 0.4:  # Check if greater than 400 mm
            pass  # Do nothing
        else:
            num_cross_members -= 1  # Reduce by one
            print("Adjusted number of cross members: Reduced by one due to free length condition.")

    elif 0.20 * cell_length_m <= free_length_old_main < 0.60 * cell_length_m:
        pass  # Do nothing

    elif 0.60 * cell_length_m <= free_length_old_main < 0.70 * cell_length_m:
        if 0.60 * cell_length_m > 0.4:
            num_cross_members += 1
            print("Adjusted number of cross members: Increased by one due to free length condition.")
        else:
           pass  # Do nothing

    elif 0.70 * cell_length_m <= free_length_old_main < 1.00 * cell_length_m:
        num_cross_members += 1
        print("Adjusted number of cross members: Increased by one due to free length condition.")

    # Recalculate the free lengths with the adjusted number of cross members
    free_length_main_member = (area_length - (num_cross_members - 1) * cell_length_m) / 2

    # Adjust number of main members based on free length of cross members
    free_length_old_cross = free_length_cross_member

    if 0 <= free_length_old_cross < 0.10 * cell_width_m:
        num_main_members -= 1
        print("Adjusted number of main members: Reduced by one due to cross member condition.")

    elif 0.10 * cell_width_m <= free_length_old_cross < 0.20 * cell_width_m:
        if 0.60 * cell_width_m > 0.4:
            pass
        else:
            num_main_members -= 1
            print("Adjusted number of main members: Reduced by one due to cross member condition.")

    elif 0.20 * cell_width_m <= free_length_old_cross < 0.60 * cell_width_m:
        pass

    elif 0.60 * cell_width_m <= free_length_old_cross < 0.70 * cell_width_m:
        if 0.60 * cell_width_m > 0.4:
            num_main_members += 1
            print("Adjusted number of main members: Increased by one due to cross member condition.")

    elif 0.70 * cell_width_m <= free_length_old_cross < 1.00 * cell_width_m:
        num_main_members += 1
        print("Adjusted number of main members: Increased by one due to cross member condition.")

    # Recalculate the free lengths with the adjusted number of main members
    free_length_cross_member = (area_width - (num_main_members - 1) * cell_width_m) / 2

    i = 0
    standard_length = 0
    while True:
        next_length = cell_length_m * (i + 1) * 1000
        if next_length > 3000:
            break
        
        i += 1
        standard_length = cell_length_m * i * 1000
    #print (" standard length= ", standard_length)
    
    n = int(area_length/3)
    main_member_df = pd.DataFrame(columns=[
        "Main Member line", "Total Length (in mm)"
        ] + [
            item for k in range(n+1) for item in [
                f"Member Name {k+1}",
                f"Length {k+1} (in mm)",
                f"First Punch {k+1} (in mm)"
            ]
        ])
    i = 1
    total_main_member = 0
    for i in range(num_main_members):
        new_row = pd.DataFrame(
            [[f"M{i + 1}", f"{area_length*1000:.2f}"]],
            columns=["Main Member line", "Total Length (in mm)"]
        )
        main_member_df = pd.concat([main_member_df, new_row], ignore_index=True)
        total_length = area_length * 1000
        solution_found = False
        tolerance = 0.01  # Acceptable rounding error  

        if total_length < 3000:                    
            main_member_df.at[main_member_df.index[-1], f"Member Name {1}"] = f"M{i + 1}_1"
            main_member_df.at[main_member_df.index[-1], f"Length {1} (in mm)"] = f"{total_length:.2f}"
            main_member_df.at[main_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(free_length_main_member - baffle_width/2):.2f}"
        else:
            A = cell_length_m * 500
            B = free_length_main_member*1000 + A

            y_max = int((3000 - A - B) // (cell_length_m * 1000))
            
            for y in reversed(range(y_max + 1)):
                first_part_length = A + B + y * cell_length_m * 1000 
                #print("First part Length", first_part_length)
                z_max = int((3000 - A - B) // (cell_length_m * 1000))
                for z in reversed(range(z_max + 1)):
                    last_part_length = A + B + z * cell_length_m * 1000
                    #print ("Last part length = ", last_part_length)  
                    remaining = total_length - first_part_length - last_part_length
                    n = remaining / standard_length
                    
                    if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                        n = round(n)
                        total_check = first_part_length + n * standard_length + last_part_length
                        if abs(total_length - total_check) <= tolerance:
                            main_member_df.at[main_member_df.index[-1], f"Member Name {1}"] = f"M{i + 1}_1"
                            main_member_df.at[main_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                            main_member_df.at[main_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{((cell_length-baffle_width)/2):.2f}"
                            
                            # Add last part (comes after n standard parts + first)
                            main_member_df.at[main_member_df.index[-1], f"Member Name {n+2}"] = f"M{i + 1}_{n+2}"
                            main_member_df.at[main_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                            main_member_df.at[main_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{((cell_length-baffle_width)/2):.2f}"      
                            for q in range(n):
                                main_member_df.at[main_member_df.index[-1], f"Member Name {q+2}"] = f"M{i + 1}_{q+2}"
                                main_member_df.at[main_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                main_member_df.at[main_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{((cell_length-baffle_width)/ 2):.2f}"                                      
                            total_main_member = total_main_member + n + 2
                            solution_found = True                                    
                            break

                if solution_found:
                    break

            if not solution_found:
                print("No valid combination found.")
                
	
	
    standard_length = 0
    i = 0
    while True:
        next_length = cell_width_m * (i + 1) * 1000
        if next_length > 3000:
            break
        
        i += 1
        standard_length = cell_width_m * i * 1000
    #print (" standard length= ", standard_length)
    
    n = int(area_width/3)
    cross_member_df = pd.DataFrame(columns=[
        "Cross Member line", "Total Length (in mm)"
        ] + [
            item for k in range(n+1) for item in [
                f"Member Name {k+1}",
                f"Length {k+1} (in mm)",
                f"First Punch {k+1} (in mm)"
            ]
        ])
    i = 1
    total_cross_member = 0
    for i in range(num_cross_members):
        new_row = pd.DataFrame(
            [[f"C{i + 1}", f"{area_width*1000:.2f}"]],
            columns=["Cross Member line", "Total Length (in mm)"]
        )
        cross_member_df = pd.concat([cross_member_df, new_row], ignore_index=True)
        total_length = area_width * 1000
        solution_found = False
        tolerance = 0.01  # Acceptable rounding error  

        if total_length < 3000:                    
            cross_member_df.at[cross_member_df.index[-1], f"Member Name {1}"] = f"C{i + 1}_1"
            cross_member_df.at[cross_member_df.index[-1], f"Length {1} (in mm)"] = f"{total_length:.2f}"
            cross_member_df.at[cross_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(free_length_cross_member - baffle_width/2):.2f}"
        else:
            A = cell_width_m * 500
            B = free_length_cross_member*1000 + A
            y_max = int((3000 - A - B) // (cell_width_m * 1000))
            
            for y in reversed(range(y_max + 1)):
                first_part_length = A + B + y * cell_width_m * 1000 
                z_max = int((3000 - A - B) // (cell_width_m * 1000))
                for z in reversed(range(z_max + 1)):
                    last_part_length = A + B + z * cell_length_m * 1000
                    #print ("Last part length = ", last_part_length)  
                    remaining = total_length - first_part_length - last_part_length
                    n = remaining / standard_length
                    
                    if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                        n = round(n)
                        total_check = first_part_length + n * standard_length + last_part_length
                        if abs(total_length - total_check) <= tolerance:
                            cross_member_df.at[cross_member_df.index[-1], f"Member Name {1}"] = f"C{i + 1}_1"
                            cross_member_df.at[cross_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                            cross_member_df.at[cross_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{((cell_width - baffle_width)/2):.2f}"
                            
                            # Add last part (comes after n standard parts + first)
                            cross_member_df.at[cross_member_df.index[-1], f"Member Name {n+2}"] = f"C{i + 1}_{n+2}"
                            cross_member_df.at[cross_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                            cross_member_df.at[cross_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_width - baffle_width/2):.2f}"      
                            for q in range(n):
                                cross_member_df.at[cross_member_df.index[-1], f"Member Name {q+2}"] = f"C{i + 1}_{q+2}"
                                cross_member_df.at[cross_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                cross_member_df.at[cross_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_width - baffle_width/ 2):.2f}"                                      
                            total_cross_member = total_cross_member + n + 2
                            solution_found = True                                    
                            break

                if solution_found:
                    break

            if not solution_found:
                print("No valid combination found.")

     # Calculate total lengths after adjustments
    total_main_length = num_main_members * area_length
    total_cross_length = num_cross_members * area_width

    # Calculate running meter per square meter
    area = area_length * area_width
    running_meter_per_sqm = (total_main_length + total_cross_length) / area

      # Generate grid details text
    quadrilaterals_open_grid_details = f"""
    **Grid Details:**
    - Total Number of Main Member Lines: {num_main_members}
    - Total Number of Cross Member Lines: {num_cross_members}
    - Free Length of Main Member: {free_length_main_member * 1000:.2f} mm
    - Free Length of Cross Member: {free_length_cross_member * 1000:.2f} mm
    - Grid Size Opening: {cell_length - baffle_width} mm × {cell_width - baffle_width} mm
    - Total Length of Main Members: {total_main_length:.2f} meters
    - Total Length of Cross Members: {total_cross_length:.2f} meters
    - Running Meter per Square Meter: {running_meter_per_sqm:.2f} running meters per square meter
    """

#def plot_grid(num_main_members, num_cross_members, area_width, area_length, free_length_cross_member, free_length_main_member, cell_width_m, cell_length_m):
    fig, ax = plt.subplots(figsize=(8.27, 11.69))
    doc = ezdxf.new()
    msp = doc.modelspace()

    for i in range(num_main_members):
        y = free_length_cross_member + i * cell_width_m
        ax.plot([0, area_length], [y, y], color='orange', label='Main Member' if i == 0 else "")
        msp.add_line((0, y*1000), (area_length*1000, y*1000    ))  # Save line to DXF

    for j in range(num_cross_members):
        x = free_length_main_member + j * cell_length_m
        ax.plot([x, x], [0, area_width], color='green', label='Cross Member' if j == 0 else "")
        msp.add_line((x*1000, 0), (x*1000, area_width*1000    ))  # Save line to DXF

    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(0, area_length + 1)
    ax.set_ylim(0, area_width + 1)
    ax.set_title(f"Grid Pattern of {cell_length} mm × {cell_width} mm ")
    ax.set_xlabel(f"Length of Area (m): {area_length} m")
    ax.set_ylabel(f"Width of Area (m): {area_width} m")

    ax.set_xticks([])  # Hide x-axis ticks
    ax.set_yticks([])  # Hide y-axis ticks
    # Hide the spines (the outer box)
    for spine in ax.spines.values():
        spine.set_visible(False)
    
    plt.grid(False)        # Disables grid lines
    legend_width = 0.5 * area_width
    legend_height = 0.5 * area_length
    # Transform data coordinates to figure coordinates
    transform = ax.transData.transform
    inv_transform = ax.transAxes.inverted().transform

    # Convert (area_length, area_width) from data coordinates to figure coordinates
    legend_x, legend_y = inv_transform(transform((area_length, area_width)))
    # Add legend with bottom-left corner at (area_length, area_width)
    ax.legend(loc='lower left', bbox_to_anchor=(legend_x, legend_y), bbox_transform=ax.transAxes, fontsize=8, frameon=True)
    plt.grid()
    doc.saveas("grid_plot.dxf")
    plt.show()

    # Define a safer directory path, e.g., within the current working directory
    output_dir = "./output"

    # Ensure the directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Set the full paths for saving
    quadrilaterals_open_pdf_path = os.path.join(output_dir, "grid_plot.pdf")
    quadrilaterals_open_image_path = os.path.join(output_dir, "grid_plot.png")
    quadrilaterals_open_dxf_path = os.path.join(output_dir, "grid_plot.dxf")

    # Save the figure as a PDF and PNG
    plt.savefig(quadrilaterals_open_pdf_path, format="pdf", bbox_inches="tight")  # Save as PDF
    plt.savefig(quadrilaterals_open_image_path, format="png", bbox_inches="tight")  # Save as PNG
    doc.saveas(quadrilaterals_open_dxf_path)  # Save DXF file
    plt.show()
    plt.close(fig)

    excel_output = BytesIO()
    with pd.ExcelWriter(excel_output, engine="openpyxl") as writer:
        main_member_df.to_excel(writer, index=False, sheet_name=f"Main Member- {total_main_member} nos")
        cross_member_df.to_excel(writer, index=False, sheet_name=f"Cross Members- {total_cross_member} nos")
        #diagonal_member_df.to_excel(writer, index=False, sheet_name=f"Diagonal Member- {total_diagonal_member} nos")
    excel_output.seek(0)

    main_member_df = main_member_df.fillna('')
    cross_member_df = cross_member_df.fillna('')
    #diagonal_member_df = diagonal_member_df.fillna('')

    # Generate HTML tables for display
    html_output = f"""
    <div style='font-family: Arial, sans-serif;'>    
        <h2>Main Members Details: {total_main_member} nos</h2>
        {main_member_df.to_html(index=False, classes='dataframe')}

        <h2>Cross Members Details: {total_cross_member} nos</h2>
        {cross_member_df.to_html(index=False, classes='dataframe')}
        
        <style>
            .dataframe {{
                width: 100%;
                border-collapse: collapse;
                margin-bottom: 20px;
                color: #000000;
                background-color: #ffffff;
            }}
            .dataframe th {{
                background-color: #ffa8B6;
                color: #373435;
                padding: 10px;
                text-align: center;
            }}
            .dataframe td {{
                padding: 8px;
                border-bottom: 1px solid #ddd;
                text-align: center;
            }}
            .dataframe tr:nth-child(even) {{
                background-color: #f2f2f2;
            }}

            /* Dark theme overrides */
            @media (prefers-color-scheme: dark) {{
                .dataframe {{
                    color: #ffffff;
                    background-color: #121212;
                }}
                .dataframe th {{
                    background-color: #444;
                    color: #ffa8b6;
                }}
                .dataframe td {{
                    border-bottom: 1px solid #444;
                }}
                .dataframe tr:nth-child(even) {{
                    background-color: #1e1e1e;
                }}
            }}
        </style>
    </div>
    """

    # Save to disk for Gradio
    excel_file_path = os.path.join(output_dir, "triangular_grid_output.xlsx")
    with open(excel_file_path, "wb") as f:
        f.write(excel_output.getvalue())
    
    # ✅ Final return
    
    return quadrilaterals_open_grid_details, quadrilaterals_open_image_path, html_output, excel_file_path, quadrilaterals_open_pdf_path, quadrilaterals_open_dxf_path

def draw_triangular_open_ended(area_length, area_width, cell_opening_size, baffle_width):
    # Convert grid dimensions from mm to meters
    cell_opening_size_m = cell_opening_size / 1000  # mm to meters
    grid_size_m = (cell_opening_size + baffle_width * math.tan(math.radians(60))) / 1000
    grid_height_m = grid_size_m * math.cos(math.radians(30))
    area_length_CtoC = area_length - (baffle_width * math.cos(math.radians(30))) / 1000
    area_width_CtoC = area_width - (baffle_width * math.cos(math.radians(60))) / 1000

    triangular_open_image_path = ""
    triangular_open_pdf_path = ""
    triangular_open_dxf_path = ""


    # Validate that grid dimensions do not exceed area dimensions
    if cell_opening_size_m >= area_width or cell_opening_size_m >= area_length:
        return "Grid dimensions must be less than or equal to area dimensions." , None, None

    # Validate that baffle width does not exceed grid dimensions
    if baffle_width >= cell_opening_size:
        return "Baffle width must be less than both grid length and width." , None, None

    # Calculate number of main members
    num_main_members = int(area_width / grid_height_m) + 1  # Horizontal lines

    # Calculate initial free lengths
    free_length_vertical_member = (area_width - (num_main_members - 1) * grid_height_m) / 2

    # Adjust number of cross members based on free length of main members
    free_length_vertical_member_old = free_length_vertical_member

    if 0 <= free_length_vertical_member_old < 0.10 * grid_height_m:
        num_main_members -= 1

    elif 0.10 * grid_height_m <= free_length_vertical_member_old < 0.20 * grid_height_m:
        if 0.60 * grid_height_m > 0.4:
            pass
        else:
            num_main_members -= 1
  
    elif 0.20 * grid_height_m <= free_length_vertical_member_old < 0.60 * grid_height_m:
        pass

    elif 0.60 * grid_height_m <= free_length_vertical_member_old < 0.70 * grid_height_m:
        if 0.60 * grid_height_m > 0.4:
            num_main_members += 1
  
    elif 0.70 * grid_height_m <= free_length_vertical_member_old < 1.00 * grid_height_m:
        num_main_members += 1

    # Recalculate the free lengths with the adjusted number of main members
    free_length_vertical_member = (area_width - (num_main_members - 1) * grid_height_m) / 2

    # Calculate total lengths after adjustments
    total_main_length = num_main_members * area_length

    #Full_cross_member_length = (area_width * 1000 - baffle_width * math.cos(math.radians(60)))/(math.cos(math.radians(30)) * 1000)
    x_f = min(area_length_CtoC, area_width_CtoC / math.tan(math.radians(60)))
    
    # Calculate y_f using the equation y_f = tan(30°) * x_f
    y_f = math.tan(math.radians(60)) * x_f

    # Compute the distance between A(0,0) and B(x_f, y_f)
    Full_cross_member_length = math.sqrt(x_f**2 + y_f**2)

    p = (area_width_CtoC - (num_main_members - 1) * grid_size_m * math.cos(math.radians(30))) / 2

  
  
    i = 0
    standard_length = 0
    while True:
        next_length = grid_size_m * (i + 1) * 1000
        if next_length > 3000:
            break
        
        i += 1
        standard_length = grid_size_m * i * 1000

    # Total Length of uncomplete cross member_a
    # First intersection will be (grid_size_m)/2 far form the center width line on top main runner
    uncomplete_cross_member_a_length = 0  # Initialize total length
    i = 0
    j = 1
        
    uncomplete_cross_member_a_num  = 0  # Initialize a uncomplete_cross_member_a_num for non-zero values
    total_cross_member = 0
    # Create an empty DataFrame to store the values
    n = int(Full_cross_member_length/3)
    # Create DataFrames for all three tables
    cross_member_df = pd.DataFrame(columns=[
        "Cross Member line", "Total Length (in mm)"
        ] + [
            item for k in range(n+1) for item in [
                f"Member Name {k+1}",
                f"Length {k+1} (in mm)",
                f"First Punch {k+1} (in mm)"
            ]
        ])
    
    c_1 = 2*grid_size_m*( 1 + i) # c_1 = 2 * (p*math.tan(math.radians(30)) + grid_size_m*0.25 - (area_length -area_length_CtoC)/2 + grid_size_m * (i))  # Initialize the first value for c
    while c_1 < Full_cross_member_length:  # Loop until c exceeds the full length
        i += 1  # Increment i
        uncomplete_cross_member_a_length += c_1  # Add c to the total length
        if c_1 > 0:  # If c is non-zero, increment the uncomplete_cross_member_a_num
            #triangular_open_grid_details += f"- C_{j}: {c_1*1000:.2f} mm \n"
            uncomplete_cross_member_a_num += 1
            new_row = pd.DataFrame([[f"C{j}", f"{c_1 * 1000:.2f}"]], columns=["Cross Member line", "Total Length (in mm)"])
            cross_member_df = pd.concat([cross_member_df, new_row], ignore_index=True)
            total_length = c_1 * 1000
            solution_found = False
            tolerance = 0.01  # Acceptable rounding error

            if total_length < 3000:                    
                cross_member_df.at[cross_member_df.index[-1], f"Member Name {1}"] = f"C{j}_1"
                cross_member_df.at[cross_member_df.index[-1], f"Length {1} (in mm)"] = f"{c_1 * 1000:.2f}"
                cross_member_df.at[cross_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(p * 1000 - baffle_width + baffle_width/4) / math.cos(math.radians(30)):.2f}"
            else:
                A = grid_size_m * 500
                B = p * 1000 / math.cos(math.radians(30))
                C = 1.5 * grid_size_m * 1000

                y_max = int((3000 - A - B) // (grid_size_m * 1000))

                for y in reversed(range(y_max + 1)):
                    first_part_length = A + B + y * grid_size_m * 1000 

                    z_max = int((3000 - (C - B)) // (grid_size_m * 1000))
                    for z in reversed(range(z_max + 1)):
                        last_part_length = C - B + z * grid_size_m * 1000
                        remaining = total_length - first_part_length - last_part_length
                        n = remaining / standard_length

                        if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                            n = round(n)
                            total_check = first_part_length + n * standard_length + last_part_length
                            if abs(total_length - total_check) <= tolerance:
                                cross_member_df.at[cross_member_df.index[-1], f"Member Name {1}"] = f"C{j}_1"
                                cross_member_df.at[cross_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                                cross_member_df.at[cross_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(cell_opening_size/2):.2f}"
                                
                                # Add last part (comes after n standard parts + first)
                                cross_member_df.at[cross_member_df.index[-1], f"Member Name {n+2}"] = f"C{j}_{n+2}"
                                cross_member_df.at[cross_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                                cross_member_df.at[cross_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_opening_size/2):.2f}"      
                                for q in range(n):
                                    cross_member_df.at[cross_member_df.index[-1], f"Member Name {q+2}"] = f"C{j}_{q+2}"
                                    cross_member_df.at[cross_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                    cross_member_df.at[cross_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_opening_size/ 2):.2f}"  

                                total_cross_member = total_cross_member + n + 2                          
                                                                        
                                solution_found = True                                    
                                break

                    if solution_found:
                        break

                if not solution_found:
                    print("No valid combination found.")

            j += 1  # Increment j here
        c_1 = 2*grid_size_m*( 1 + i)  # Update c for the next iteration

 
    Full_cross_member_num = max(int(((area_length - baffle_width * (math.cos(math.radians(30))) / 1000) - uncomplete_cross_member_a_num * grid_size_m) / grid_size_m),int(area_width_CtoC /(2 * grid_height_m)) - uncomplete_cross_member_a_num)
    Full_cross_member_total_length = Full_cross_member_length * Full_cross_member_num
    Right_free_length = (area_length_CtoC - (grid_size_m - p / math.tan(math.radians(60)))) - (int((area_length_CtoC - (grid_size_m - p / math.tan(math.radians(60))))/grid_size_m))*grid_size_m
    #print (Right_free_length)

    if area_length >= area_width / (3 ** 0.5): 
        for i in range(Full_cross_member_num):
            new_row = pd.DataFrame(
                [[f"C{i + 1 + uncomplete_cross_member_a_num}", f"{Full_cross_member_length * 1000:.2f}"]],
                columns=["Cross Member line", "Total Length (in mm)"]
            )
            cross_member_df = pd.concat([cross_member_df, new_row], ignore_index=True)
            total_length = Full_cross_member_length * 1000
            solution_found = False
            tolerance = 0.01  # Acceptable rounding error  

            if total_length < 3000:                    
                cross_member_df.at[cross_member_df.index[-1], f"Member Name {1}"] = f"C{i + 1 + uncomplete_cross_member_a_num}_1"
                cross_member_df.at[cross_member_df.index[-1], f"Length {1} (in mm)"] = f"{total_length:.2f}"
                cross_member_df.at[cross_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(p * 1000 - baffle_width + baffle_width/4) / math.cos(math.radians(30)):.2f}"
            else:
                A = grid_size_m * 500
                B = p * 1000 / math.cos(math.radians(30))
                C = 1.5 * grid_size_m * 1000

                y_max = int((3000 - A - B) // (grid_size_m * 1000))
                
                for y in reversed(range(y_max + 1)):
                    first_part_length = A + B + y * grid_size_m * 1000                        

                    z_max = int((3000 - A - B) // (grid_size_m * 1000))
                    for z in reversed(range(z_max + 1)):
                        last_part_length = A + B + z * grid_size_m * 1000
                        remaining = total_length - first_part_length - last_part_length
                        n = remaining / standard_length
                        
                        if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                            n = round(n)
                            total_check = first_part_length + n * standard_length + last_part_length
                            if abs(total_length - total_check) <= tolerance:
                                cross_member_df.at[cross_member_df.index[-1], f"Member Name {1}"] = f"C{i + 1 + uncomplete_cross_member_a_num}_1"
                                cross_member_df.at[cross_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                                cross_member_df.at[cross_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(cell_opening_size/2):.2f}"
                                
                                # Add last part (comes after n standard parts + first)
                                cross_member_df.at[cross_member_df.index[-1], f"Member Name {n+2}"] = f"C{i + 1 + uncomplete_cross_member_a_num}_{n+2}"
                                cross_member_df.at[cross_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                                cross_member_df.at[cross_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_opening_size/2):.2f}"      
                                for q in range(n):
                                    cross_member_df.at[cross_member_df.index[-1], f"Member Name {q+2}"] = f"C{i + 1 + uncomplete_cross_member_a_num}_{q+2}"
                                    cross_member_df.at[cross_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                    cross_member_df.at[cross_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_opening_size/ 2):.2f}"  
                                total_cross_member = total_cross_member + n + 2                                    
                                solution_found = True                                    
                                break

                    if solution_found:
                        break

                if not solution_found:
                    print("No valid combination found.")
    else:
        for i in range(Full_cross_member_num):
            new_row = pd.DataFrame([[f"C{i + 1 + uncomplete_cross_member_a_num}", f"{Full_cross_member_length * 1000:.2f}"]], columns=["Cross Member line", "Total Length (in mm)"])
          
        
            cross_member_df = pd.concat([cross_member_df, new_row], ignore_index=True)
            total_length = Full_cross_member_length * 1000
            solution_found = False
            tolerance = 0.01  # Acceptable rounding error  

            if total_length < 3000:                    
                cross_member_df.at[cross_member_df.index[-1], f"Member Name {1}"] = f"C{i + 1 + uncomplete_cross_member_a_num}_1"
                cross_member_df.at[cross_member_df.index[-1], f"Length {1} (in mm)"] = f"{total_length:.2f}"
                cross_member_df.at[cross_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(p * 1000 - baffle_width + baffle_width/4) / math.cos(math.radians(30)):.2f}"
            else:
                A = grid_size_m * 500
                B = p * 1000 / math.cos(math.radians(30))
                C = 1.5 * grid_size_m * 1000


                y_max = int((3000 - A - 2*Right_free_length*1000) // (grid_size_m * 1000))

                for y in reversed(range(y_max + 1)):
                    first_part_length = 2*Right_free_length*1000 + A + y * grid_size_m * 1000                        

                    z_max = int((3000 - (C - B)) // (grid_size_m * 1000))
                    for z in reversed(range(z_max + 1)):
                        last_part_length = C - B + z * grid_size_m * 1000
                        remaining = total_length - first_part_length - last_part_length
                        n = remaining / standard_length

                        if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                            n = round(n)
                            total_check = first_part_length + n * standard_length + last_part_length
                            if abs(total_length - total_check) <= tolerance:
                                cross_member_df.at[cross_member_df.index[-1], f"Member Name {1}"] = f"C{i + 1 + uncomplete_cross_member_a_num}_1"
                                cross_member_df.at[cross_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                                cross_member_df.at[cross_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(cell_opening_size/2):.2f}"
                                
                                # Add last part (comes after n standard parts + first)
                                cross_member_df.at[cross_member_df.index[-1], f"Member Name {n+2}"] = f"C{i + 1 + uncomplete_cross_member_a_num}_{n+2}"
                                cross_member_df.at[cross_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                                cross_member_df.at[cross_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_opening_size/2):.2f}"      
                                for q in range(n):
                                    cross_member_df.at[cross_member_df.index[-1], f"Member Name {q+2}"] = f"C{i + 1 + uncomplete_cross_member_a_num}_{q+2}"
                                    cross_member_df.at[cross_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                    cross_member_df.at[cross_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_opening_size/ 2):.2f}"                                                               
                                total_cross_member = total_cross_member + n + 2
                                solution_found = True                                    
                                break

                    if solution_found:
                        break

                if not solution_found:
                    print("No valid combination found.")


    # Total Length of uncomplete cross member_b Bottom right
    uncomplete_cross_member_b_length = 0  # Initialize total length
    i = 0
    j = 1
    # triangular_open_grid_details += f"""
    # Uncomplete Cross Member (Bottom Right) Details:
    # """
    slope =  math.tan(math.radians(60))
    uncomplete_cross_member_b_num = 0  # Initialize a uncomplete_cross_member_b_num for non-zero values
    c_2 = ((area_length_CtoC - (((((area_width - area_width_CtoC)/2) - (area_width - (area_width -area_width_CtoC)/2))/slope + (area_length - area_length_CtoC)/2 + (uncomplete_cross_member_a_num + 1) * grid_size_m + (Full_cross_member_num - 1)*grid_size_m) - (area_length - area_length_CtoC)/2)) - (i+1) * grid_size_m)*2
    # (area_length_CtoC - ((uncomplete_cross_member_a_num + Full_cross_member_num) * grid_size_m + 0.5 * abs(grid_size_m - 2 * p / math.cos(math.radians(30)))) - grid_size_m + i * grid_size_m) * 2  # Initialize the first value for c  # Initialize the first value for c
      # (area_length_CtoC - ((uncomplete_cross_member_a_num + Full_cross_member_num) * grid_size_m + (grid_size_m - 2 * p * math.tan(math.radians(30)))) - grid_size_m + i * grid_size_m) * 2                #

    while c_2 > (area_length - area_length_CtoC)/2:  # Loop until c exceeds the full length
        i += 1  # Increment i

        if 0 < c_2 < Full_cross_member_length:  # If c is non-zero, increment the uncomplete_cross_member_b_num
            #triangular_open_grid_details += f" - C_{j+Full_cross_member_num+uncomplete_cross_member_a_num}: {c_2*1000:.2f} mm \n"
            uncomplete_cross_member_b_length += c_2  # Add c to the total length
            uncomplete_cross_member_b_num += 1  # Increment the counter
            new_row = pd.DataFrame([[f"C{j + Full_cross_member_num + uncomplete_cross_member_a_num}", f"{c_2 * 1000:.2f}"]], columns=["Cross Member line", "Total Length (in mm)"])
            cross_member_df = pd.concat([cross_member_df, new_row], ignore_index=True)
            total_length = c_2 * 1000
            solution_found = False
            tolerance = 0.01  # Acceptable rounding error

            if total_length < 3000:                    
                cross_member_df.at[cross_member_df.index[-1], f"Member Name {1}"] = f"C{j + Full_cross_member_num + uncomplete_cross_member_a_num}_1"
                cross_member_df.at[cross_member_df.index[-1], f"Length {1} (in mm)"] = f"{c_2 * 1000:.2f}"
                cross_member_df.at[cross_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(p * 1000 - baffle_width + baffle_width/4) / math.cos(math.radians(30)):.2f}"
            else:
                A = grid_size_m * 500
                B = p * 1000 / math.cos(math.radians(30))
                C = 1.5 * grid_size_m * 1000    

                y_max = int((3000 - A - Right_free_length*2000) // (grid_size_m * 1000))

                for y in reversed(range(y_max + 1)):
                    first_part_length = A + Right_free_length*2000 + y * grid_size_m * 1000                        

                    z_max = int((3000 - (A + B)) // (grid_size_m * 1000))
                    for z in reversed(range(z_max + 1)):
                        last_part_length = A +  B + z * grid_size_m * 1000
                        remaining = total_length - first_part_length - last_part_length
                        n = remaining / standard_length

                        if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                            n = round(n)
                            total_check = first_part_length + n * standard_length + last_part_length
                            if abs(total_length - total_check) <= tolerance:
                                cross_member_df.at[cross_member_df.index[-1], f"Member Name {1}"] = f"C{j + Full_cross_member_num + uncomplete_cross_member_a_num}_1"
                                cross_member_df.at[cross_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                                cross_member_df.at[cross_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(cell_opening_size/2):.2f}"
                                
                                # Add last part (comes after n standard parts + first)
                                cross_member_df.at[cross_member_df.index[-1], f"Member Name {n+2}"] = f"C{j + Full_cross_member_num + uncomplete_cross_member_a_num}_{n+2}"
                                cross_member_df.at[cross_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                                cross_member_df.at[cross_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_opening_size/2):.2f}"      
                                for q in range(n):
                                    cross_member_df.at[cross_member_df.index[-1], f"Member Name {q+2}"] = f"C{j + Full_cross_member_num + uncomplete_cross_member_a_num}_{q+2}"
                                    cross_member_df.at[cross_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                    cross_member_df.at[cross_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_opening_size/ 2):.2f}"                            
                                total_cross_member = total_cross_member + n + 2                                        
                                solution_found = True                                    
                                break

                    if solution_found:
                        break
                if not solution_found:
                    print("No valid combination found.")
            j += 1  # Increment j here
        c_2 =  ((area_length_CtoC - (((((area_width - area_width_CtoC)/2) - (area_width - (area_width -area_width_CtoC)/2))/slope + (area_length - area_length_CtoC)/2 + (uncomplete_cross_member_a_num + 1) * grid_size_m + (Full_cross_member_num - 1)*grid_size_m) - (area_length - area_length_CtoC)/2)) - (i+1) * grid_size_m)*2  # Update c for the next iteration

    # Print the number of non-zero iterations and the total length after the loop
    # triangular_open_grid_details += f"""    
    # - Number of uncomplete cross members in Bottom Right side: {uncomplete_cross_member_b_num}
    # - Total Length of uncomplete cross member in Bottom Right side: {uncomplete_cross_member_b_length:.2f} meters
    # """
    # ===== TABLE 2: Diagonal Member line =====
    n = int(Full_cross_member_length/3)
    diagonal_member_df = pd.DataFrame(columns=[
        "Diagonal Member line", "Total Length (in mm)"
        ] + [
            item for k in range(n+1) for item in [
                f"Member Name {k+1}",
                f"Length {k+1} (in mm)",
                f"First Punch {k+1} (in mm)"
            ]
        ])
    
    Full_diagonal_member_length = Full_cross_member_length #(area_width * 1000 - baffle_width * math.cos(math.radians(60))) / (math.cos(math.radians(30)) * 1000)
    total_diagonal_member = 0

    # Total Length of uncomplete diagonal member_a Top Left
    uncomplete_diagonal_member_a_length = 0  # Initialize total length
    i = 0
    j = 1
    # triangular_open_grid_details += f"""
    # Uncomplete Diagonal Member (Top Right) Details:
    # """    
    uncomplete_diagonal_member_a_num = 0  # Initialize a uncomplete_diagonal_member_a_num for non-zero values
    c_3 =  2*(area_length_CtoC-(int((area_length_CtoC-(grid_size_m-2*p/math.tan(math.radians(60))))/grid_size_m))*grid_size_m-(grid_size_m-2*p/math.tan(math.radians(60))) + i*grid_size_m)  # Initialize the first value for c
    while c_3 < Full_diagonal_member_length:  # Loop until c exceeds the full length
        i += 1  # Increment i
        if c_3 > 0:  # If c is non-zero, increment the uncomplete_diagonal_member_a_num
            #triangular_open_grid_details += f"- D_{j}: {c_3*1000:.2f} mm \n"
            new_row = pd.DataFrame([[f"D{j}", f"{c_3 * 1000:.2f}"]], columns=["Diagonal Member line", "Total Length (in mm)"])
            diagonal_member_df = pd.concat([diagonal_member_df, new_row], ignore_index=True)
            total_length = c_3 * 1000
            solution_found = False
            tolerance = 0.01  # Acceptable rounding error

            if total_length < 3000:                    
                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {1}"] = f"D{j}_1"
                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {1} (in mm)"] = f"{c_3 * 1000:.2f}"
                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(p * 1000 - baffle_width + baffle_width/4) / math.cos(math.radians(30)):.2f}"
            else:
                A = grid_size_m * 500
                B = p * 1000 / math.cos(math.radians(30))
                C = 1.5 * grid_size_m * 1000

                y_max = int((3000 - A - B) // (grid_size_m * 1000))
                for y in reversed(range(y_max + 1)):
                    first_part_length = A + B + y * grid_size_m * 1000

                    z_max = int((3000 - A - Right_free_length*2000) // (grid_size_m * 1000))
                    for z in reversed(range(z_max + 1)):
                        last_part_length = A + Right_free_length*2000 + z * grid_size_m * 1000
                        
                        remaining = total_length - first_part_length - last_part_length
                        n = remaining / standard_length

                        if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                            n = round(n)
                            total_check = first_part_length + n * standard_length + last_part_length
                            if abs(total_length - total_check) <= tolerance:
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {1}"] = f"D{j}_1"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(cell_opening_size/2):.2f}"
                                
                                # Add last part (comes after n standard parts + first)
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {n+2}"] = f"D{j}_{n+2}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_opening_size/2):.2f}"      
                                for q in range(n):
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {q+2}"] = f"D{j}_{q+2}"
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_opening_size/ 2):.2f}"                            
                                total_diagonal_member = total_diagonal_member + n + 2                                        
                                solution_found = True                                    
                                break

                    if solution_found:
                        break

                if not solution_found:
                    print("No valid combination found.")

            uncomplete_diagonal_member_a_num += 1
            uncomplete_diagonal_member_a_length += c_3  # Add c to the total length
            j += 1
        c_3 = 2*(area_length_CtoC-(int((area_length_CtoC-(grid_size_m-2*p/math.tan(math.radians(60))))/grid_size_m))*grid_size_m-(grid_size_m-2*p/math.tan(math.radians(60))) + i*grid_size_m)  # Update c for the next iteration

    # Print the number of non-zero iterations and the total length after the loop
    # triangular_open_grid_details += f"""  
    # - Number of uncomplete diagonal members in top right side: {uncomplete_diagonal_member_a_num}
    # - Total Length of uncomplete diagonal member in top right side: {uncomplete_diagonal_member_a_length:.2f}
    # """

    #Full_diagonal_member_num = max(int(((area_length - baffle_width * (math.cos(math.radians(30))) / 1000) - uncomplete_diagonal_member_a_num * grid_size_m -(area_length_CtoC-(int((area_length_CtoC-(grid_size_m-2*p/math.tan(math.radians(60))))/grid_size_m))*grid_size_m-(grid_size_m-2*p/math.tan(math.radians(60))))) / grid_size_m) + 1, 
                                   #int((area_width_CtoC - (area_length_CtoC-(int((area_length_CtoC-(grid_size_m-2*p/math.tan(math.radians(60))))/grid_size_m))*grid_size_m-(grid_size_m-2*p/math.tan(math.radians(60))))) /(2 * grid_height_m)) - uncomplete_cross_member_a_num + 1) 


    Full_diagonal_member_num = 0
    
    if area_length >= area_width / (3 ** 0.5):
        slope = math.tan(math.radians(120))
        y_start = area_width - (area_width -area_width_CtoC)/2
        y_end = (area_width - area_width_CtoC)/2 #slope*(x_start - x_end) + y_end
        #x_end = (area_length - area_length_CtoC)/2 + (uncomplete_cross_member_a_num + 1) * grid_size_m       #(p*math.tan(math.radians(30)) + grid_size_m*0.25) + (area_length -area_length_CtoC)/2 + (uncomplete_cross_member_a_num) * grid_size_m
        x_start = (area_length - area_length_CtoC)/2 + grid_size_m - 2*p*math.tan(math.radians(30)) #(y_start - y_end)/slope + x_end
        x_end = (y_end - y_start)/slope + x_start
        first_iteration = True
 
        while x_end < area_length - (area_length - area_length_CtoC)/2 :
            
            # Plot the line with the current x_start and y_start values
            # ax.plot([x_start, x_end], [y_start, y_end], color='cyan') #, label='Full diagonal Member' if first_iteration else "")
            # msp.add_line((x_start*1000, y_start*1000), (x_end*1000, y_end*1000))  # Save line to DXF
            Full_diagonal_member_num = Full_diagonal_member_num + 1
            x_start = x_start + grid_size_m
            x_end = (y_end - y_start)/slope + x_start
            
    else:
        slope = math.tan(math.radians(120))
        x_start = (area_length - area_length_CtoC)/2
        x_end = (area_length + area_length_CtoC)/2
        y_start = (num_main_members-1) * grid_height_m + (area_width -area_width_CtoC)/2
        y_end = y_start + (x_end-x_start)*slope

        while y_end >= (area_width -area_width_CtoC)/2:
            Full_diagonal_member_num = Full_diagonal_member_num + 1
            y_start = y_start - 2 * grid_height_m
            y_end = y_start + (x_end-x_start)*slope        


            
    #(int((area_length_CtoC-(grid_size_m-2*p/math.tan(math.radians(60))))/grid_size_m)) - uncomplete_diagonal_member_a_num    
    Full_diagonal_member_total_length = Full_diagonal_member_length * Full_diagonal_member_num

    if area_length >= area_width / (3 ** 0.5): 
        for i in range(Full_diagonal_member_num):
            new_row = pd.DataFrame(
                [[f"D{i + 1 + uncomplete_diagonal_member_a_num}", f"{Full_diagonal_member_length * 1000:.2f}"]],
                columns=["Diagonal Member line", "Total Length (in mm)"]
            )
            diagonal_member_df = pd.concat([diagonal_member_df, new_row], ignore_index=True)
            total_length = Full_diagonal_member_length * 1000
            solution_found = False
            tolerance = 0.01  # Acceptable rounding error  

            if total_length < 3000:                    
                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {1}"] = f"D{i + 1 + uncomplete_diagonal_member_a_num}_1"
                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {1} (in mm)"] = f"{total_length:.2f}"
                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(p * 1000 - baffle_width + baffle_width/4) / math.cos(math.radians(30)):.2f}"
            else:
                A = grid_size_m * 500
                B = p * 1000 / math.cos(math.radians(30))
                C = 1.5 * grid_size_m * 1000

                y_max = int((3000 - A - B) // (grid_size_m * 1000))
                
                for y in reversed(range(y_max + 1)):
                    first_part_length = A + B + y * grid_size_m * 1000                        

                    z_max = int((3000 - A - B) // (grid_size_m * 1000))
                    for z in reversed(range(z_max + 1)):
                        last_part_length = A + B + z * grid_size_m * 1000
                        remaining = total_length - first_part_length - last_part_length
                        n = remaining / standard_length
                        
                        if 0 <= last_part_length < 3000 and n >= 0 and abs(n - round(n)) <= 1e-6:
                            n = round(n)
                            total_check = first_part_length + n * standard_length + last_part_length
                            if abs(total_length - total_check) <= tolerance:
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {1}"] = f"D{i + 1 + uncomplete_diagonal_member_a_num}_1"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(cell_opening_size/2):.2f}"
                                
                                # Add last part (comes after n standard parts + first)
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {n+2}"] = f"D{i + 1 + uncomplete_diagonal_member_a_num}_{n+2}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_opening_size/2):.2f}"      
                                for q in range(n):
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {q+2}"] = f"D{i + 1 + uncomplete_diagonal_member_a_num}_{q+2}"
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_opening_size/ 2):.2f}"                                      
                                total_diagonal_member = total_diagonal_member + n + 2
                                solution_found = True                                    
                                break

                    if solution_found:
                        break

                if not solution_found:
                    print("No valid combination found.")
    else:
        for i in range(Full_diagonal_member_num):
            new_row = pd.DataFrame([[f"C{i + 1 + uncomplete_diagonal_member_a_num}", f"{Full_diagonal_member_length * 1000:.2f}"]], columns=["Diagonal Member line", "Total Length (in mm)"])
          
        
            diagonal_member_df = pd.concat([diagonal_member_df, new_row], ignore_index=True)
            total_length = Full_diagonal_member_length * 1000
            solution_found = False
            tolerance = 0.01  # Acceptable rounding error  

            if total_length < 3000:                    
                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {1}"] = f"D{i + 1 + uncomplete_diagonal_member_a_num}_1"
                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {1} (in mm)"] = f"{total_length:.2f}"
                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(p * 1000 - baffle_width + baffle_width/4) / math.cos(math.radians(30)):.2f}"
            else:
                A = grid_size_m * 500
                B = p * 1000 / math.cos(math.radians(30))
                C = 1.5 * grid_size_m * 1000


                y_max = int((3000 - (C - B)) // (grid_size_m * 1000))

                for y in reversed(range(y_max + 1)):
                    first_part_length = C - B + y * grid_size_m * 1000                        

                    z_max = int((3000 - A - 2*Right_free_length*1000) // (grid_size_m * 1000))
                    for z in reversed(range(z_max + 1)):
                        last_part_length = 2*Right_free_length*1000 + A + z * grid_size_m * 1000
                        remaining = total_length - first_part_length - last_part_length
                        n = remaining / standard_length

                        if 0 <= last_part_length < 3000 and n >= 0 and abs(n - round(n)) <= 1e-6:
                            n = round(n)
                            total_check = first_part_length + n * standard_length + last_part_length
                            if abs(total_length - total_check) <= tolerance:
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {1}"] = f"D{i + 1 + uncomplete_diagonal_member_a_num}_1"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(cell_opening_size/2):.2f}"
                                
                                # Add last part (comes after n standard parts + first)
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {n+2}"] = f"D{i + 1 + uncomplete_diagonal_member_a_num}_{n+2}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_opening_size/2):.2f}"      
                                for q in range(n):
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {q+2}"] = f"D{i + 1 + uncomplete_diagonal_member_a_num}_{q+2}"
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_opening_size/ 2):.2f}"                                                               
                                total_diagonal_member = total_diagonal_member + n + 2
                                solution_found = True                                    
                                break

                    if solution_found:
                        break

                if not solution_found:
                    print("No valid combination found.")





    # Total Length of uncomplete diagonal member_b Left Bottom
    uncomplete_diagonal_member_b_length = 0  # Initialize total length
    i = 0
    j = 1
    # triangular_open_grid_details += f"""
    # Uncomplete Diagonal Member (Bottom left) Details:
    # """ 
    slope = - math.tan(math.radians(60))
    uncomplete_diagonal_member_b_num = 0  # Initialize a uncomplete_diagonal_member_b_num for non-zero values

    c_4 = (((((area_width - area_width_CtoC)/2) - (area_width - (area_width -area_width_CtoC)/2))/slope + ((area_length - area_length_CtoC)/2 + grid_size_m - 2*p*math.tan(math.radians(30))) - grid_size_m) - (area_length - area_length_CtoC)/2 - i * grid_size_m)*2
     #((grid_size_m-2*p/math.tan(math.radians(60)))+abs(grid_size_m - 2*p/math.cos(math.radians(30)))*0.5 + i * grid_size_m)*2 # Initialize the first value for
     #(area_length_CtoC - ((uncomplete_cross_member_a_num + Full_cross_member_num) * grid_size_m + (grid_size_m - 2 * p * math.tan(math.radians(30)))) - grid_size_m + i * grid_size_m) * 2
    while c_4 > (area_length - area_length_CtoC)/2:  # Loop until c exceeds the full length
        i += 1  # Increment i

        if 0 < c_4 < Full_diagonal_member_length:  # If c_4 is non-zero, increment the uncomplete_diagonal_member_b_num
            #triangular_open_grid_details += f"- D_{j+Full_diagonal_member_num+uncomplete_diagonal_member_a_num}: {c_4*1000:.2f} mm\n"
            uncomplete_diagonal_member_b_length += c_4  # Add c to the total length
            uncomplete_diagonal_member_b_num += 1  # Increment the counter
            new_row = pd.DataFrame([[f"D{j + Full_diagonal_member_num + uncomplete_diagonal_member_a_num}", f"{c_4 * 1000:.2f}"]], columns=["Diagonal Member line", "Total Length (in mm)"])
            diagonal_member_df = pd.concat([diagonal_member_df, new_row], ignore_index=True)
            total_length = c_4 * 1000
            solution_found = False
            tolerance = 0.01  # Acceptable rounding error

            if total_length < 3000:                    
                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {1}"] = f"D{j + Full_diagonal_member_num + uncomplete_diagonal_member_a_num}_1"
                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {1} (in mm)"] = f"{c_4 * 1000:.2f}"
                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(p * 1000 - baffle_width + baffle_width/4) / math.cos(math.radians(30)):.2f}"
            else:
                A = grid_size_m * 500
                B = p * 1000 / math.cos(math.radians(30))
                C = 1.5 * grid_size_m * 1000    

                y_max = int((3000 - C + B) // (grid_size_m * 1000))

                for y in reversed(range(y_max + 1)):
                    first_part_length = C - B + y * grid_size_m * 1000                        

                    z_max = int((3000 - (A + B)) // (grid_size_m * 1000))
                    for z in reversed(range(z_max + 1)):
                        last_part_length = A +  B + z * grid_size_m * 1000
                        remaining = total_length - first_part_length - last_part_length
                        n = remaining / standard_length

                        if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                            n = round(n)
                            total_check = first_part_length + n * standard_length + last_part_length
                            if abs(total_length - total_check) <= tolerance:
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {1}"] = f"D{j + Full_diagonal_member_num + uncomplete_diagonal_member_a_num}_1"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(cell_opening_size/2):.2f}"
                                
                                # Add last part (comes after n standard parts + first)
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {n+2}"] = f"D{j + Full_diagonal_member_num + uncomplete_diagonal_member_a_num}_{n+2}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_opening_size/2):.2f}"      
                                for q in range(n):
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {q+2}"] = f"D{j + Full_diagonal_member_num + uncomplete_diagonal_member_a_num}_{q+2}"
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_opening_size/ 2):.2f}"                            
                                                                        
                                total_diagonal_member = total_diagonal_member + n + 2
                                solution_found = True                                    
                                break

                    if solution_found:
                        break
                if not solution_found:
                    print("No valid combination found.")
            j += 1

        c_4 = (((((area_width - area_width_CtoC)/2) - (area_width - (area_width -area_width_CtoC)/2))/slope + ((area_length - area_length_CtoC)/2 + grid_size_m - 2*p*math.tan(math.radians(30))) - grid_size_m) - (area_length - area_length_CtoC)/2 - i * grid_size_m)*2  # Update c for the next iteration
        # (area_length_CtoC - ((uncomplete_cross_member_a_num + Full_cross_member_num) * grid_size_m + (grid_size_m - 2 * p * math.tan(math.radians(30)))) - grid_size_m + i * grid_size_m) * 2  #

    # Print the number of non-zero iterations and the total length after the loop
    # triangular_open_grid_details += f"""    
    # - Number of uncomplete diagonal members in Bottom left side: {uncomplete_diagonal_member_b_num}
    # - Total Length of uncomplete diagonal member in Bottom left side: {uncomplete_diagonal_member_b_length:.2f}
    # """

    Total_running_meter = total_main_length + Full_cross_member_total_length + uncomplete_cross_member_a_length + uncomplete_cross_member_b_length + Full_diagonal_member_total_length + uncomplete_diagonal_member_a_length + uncomplete_diagonal_member_b_length

   
    # ===== TABLE 3:Main Members Details ===== 
    n = int(area_length/3)
    main_member_df = pd.DataFrame(columns=[
        "Main Member line", "Total Length (in mm)"
        ] + [
            item for k in range(n+1) for item in [
                f"Member Name {k+1}",
                f"Length {k+1} (in mm)",
                f"First Punch {k+1} (in mm)"
            ]
        ])
    i = 1
    total_main_member = 0
    for i in range(num_main_members):
        if i % 2 == 0:
            new_row = pd.DataFrame(
                [[f"M{i + 1}", f"{area_length*1000:.2f}"]],
                columns=["Main Member line", "Total Length (in mm)"]
            )
            main_member_df = pd.concat([main_member_df, new_row], ignore_index=True)
            total_length = area_length * 1000
            solution_found = False
            tolerance = 0.01  # Acceptable rounding error  

            if total_length < 3000:                    
                main_member_df.at[main_member_df.index[-1], f"Member Name {1}"] = f"M{i + 1}_1"
                main_member_df.at[main_member_df.index[-1], f"Length {1} (in mm)"] = f"{total_length:.2f}"
                main_member_df.at[main_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(grid_size_m*1000 - (p * 1000 + 1.5 * baffle_width) * math.tan(math.radians(30)) + (area_length - area_length_CtoC)/2):.2f}"
            else:
                A = grid_size_m * 500
                B = p * 1000 / math.cos(math.radians(30))
                C = 1.5 * grid_size_m * 1000
                D = grid_size_m * 1000 - p * 1000 * math.tan(math.radians(30)) + (area_length - area_length_CtoC)*1000/2
                E = Right_free_length*1000 + (area_length - area_length_CtoC)*1000/2

                y_max = int((3000 - A - D) // (grid_size_m * 1000))
                
                for y in reversed(range(y_max + 1)):
                    first_part_length = A + D + y * grid_size_m * 1000 
                    z_max = int((3000 - A - E) // (grid_size_m * 1000))
                    for z in reversed(range(z_max + 1)):
                        last_part_length = A + E + z * grid_size_m * 1000
                        #print ("Last part length = ", last_part_length)  
                        remaining = total_length - first_part_length - last_part_length
                        n = remaining / standard_length
                        
                        if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                            n = round(n)
                            total_check = first_part_length + n * standard_length + last_part_length
                            if abs(total_length - total_check) <= tolerance:
                                main_member_df.at[main_member_df.index[-1], f"Member Name {1}"] = f"M{i + 1}_1"
                                main_member_df.at[main_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                                main_member_df.at[main_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(cell_opening_size/2):.2f}"
                                
                                # Add last part (comes after n standard parts + first)
                                main_member_df.at[main_member_df.index[-1], f"Member Name {n+2}"] = f"M{i + 1}_{n+2}"
                                main_member_df.at[main_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                                main_member_df.at[main_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_opening_size/2):.2f}"      
                                for q in range(n):
                                    main_member_df.at[main_member_df.index[-1], f"Member Name {q+2}"] = f"M{i + 1}_{q+2}"
                                    main_member_df.at[main_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                    main_member_df.at[main_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_opening_size/ 2):.2f}"                                      
                                total_main_member = total_main_member + n + 2
                                solution_found = True                                    
                                break

                    if solution_found:
                        break

                if not solution_found:
                    print("No valid combination found.")
        else:
            new_row = pd.DataFrame(
                [[f"M{i + 1}", f"{area_length*1000:.2f}"]],
                columns=["Main Member line", "Total Length (in mm)"]
            )
            main_member_df = pd.concat([main_member_df, new_row], ignore_index=True)
            total_length = area_length * 1000
            solution_found = False
            tolerance = 0.01  # Acceptable rounding error  

            if total_length < 3000:                    
                main_member_df.at[main_member_df.index[-1], f"Member Name {1}"] = f"M{i + 1}_1"
                main_member_df.at[main_member_df.index[-1], f"Length {1} (in mm)"] = f"{total_length:.2f}"
                main_member_df.at[main_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(grid_size_m*1000 - (p * 1000 + 1.5 * baffle_width) * math.tan(math.radians(30)) + (area_length - area_length_CtoC)/2):.2f}"
            else:
                A = grid_size_m * 500
                B = p * 1000 / math.cos(math.radians(30))
                C = 1.5 * grid_size_m * 1000
                D = grid_size_m * 1000 - p * 1000 * math.tan(math.radians(30)) + (area_length - area_length_CtoC)*1000/2 - grid_size_m * 500
                E = Right_free_length*1000 + (area_length - area_length_CtoC)*1000/2 + grid_size_m * 500

                y_max = int((3000 - A - D) // (grid_size_m * 1000))
                
                for y in reversed(range(y_max + 1)):
                    first_part_length = A + D + y * grid_size_m * 1000 

                    z_max = int((3000 - A - E) // (grid_size_m * 1000))
                    for z in reversed(range(z_max + 1)):
                        last_part_length = A + E + z * grid_size_m * 1000  
                        remaining = total_length - first_part_length - last_part_length
                        n = remaining / standard_length
                        
                        if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                            n = round(n)
                            total_check = first_part_length + n * standard_length + last_part_length
                            if abs(total_length - total_check) <= tolerance:
                                main_member_df.at[main_member_df.index[-1], f"Member Name {1}"] = f"M{i + 1}_1"
                                main_member_df.at[main_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                                main_member_df.at[main_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(cell_opening_size/2):.2f}"
                                
                                # Add last part (comes after n standard parts + first)
                                main_member_df.at[main_member_df.index[-1], f"Member Name {n+2}"] = f"M{i + 1}_{n+2}"
                                main_member_df.at[main_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                                main_member_df.at[main_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_opening_size/2):.2f}"      
                                for q in range(n):
                                    main_member_df.at[main_member_df.index[-1], f"Member Name {q+2}"] = f"M{i + 1}_{q+2}"
                                    main_member_df.at[main_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                    main_member_df.at[main_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_opening_size/ 2):.2f}"                                      
                                total_main_member = total_main_member + n + 2
                                solution_found = True                                    
                                break

                    if solution_found:
                        break

                if not solution_found:
                    print("No valid combination found.")
    
    # Initialize output string
    triangular_open_grid_details = "**Grid Details:**\n"
    # Print results
    
    # - Total Number of Main Member Lines: {num_main_members}
    # - Length of Main Member Lines: {area_length*1000:.2f} mm = {area_length:.2f} meters
    # - Total Length of Main Member Lines: {total_main_length:.2f} meters 
    
    # - Total Number of Full Cross Member Lines: {Full_cross_member_num} 
    # - Length of Full Cross Member Lines: {Full_cross_member_length*1000:.2f} mm = {Full_cross_member_length:.2f} meters 
    # - Total Length of Full Cross Member Lines: {Full_cross_member_total_length:.2f} meters
    
    # - Total Number of Full diagonal Member Lines: {Full_diagonal_member_num} 
    # - Length of Full diagonal Member Lines: {Full_diagonal_member_length*1000:.2f} mm = {Full_diagonal_member_length:.2f} meters 
    # - Total Length of Full diagonal Member Lines: {Full_diagonal_member_total_length:.2f} meters 
   
    triangular_open_grid_details += f"""
    - triangular Grid Size Opening (in mm): {cell_opening_size} mm × {cell_opening_size} mm × {cell_opening_size} mm
    - Total required Running meter: {Total_running_meter:.2f} meters
    - Running Meter per Square Meter: {Total_running_meter/(area_length*area_width):.2f} running meters per square meter 
    
    - Main member Standard Punching Dimensions: {(cell_opening_size - 0.4):.2f} mm
    - Diaginal member Standard Punching Dimensions: {(cell_opening_size - 0.4):.2f} mm
    - Cross member Standard Punching Dimensions: {(grid_size_m*1000 - 29.5):.2f} mm 
    """

    

    # Return all necessary values
    # return ( grid_height_m, grid_size_m, uncomplete_cross_member_a_length, area_length_CtoC, p, area_width_CtoC)

  #def plot_grid(grid_size_m, num_main_members, total_main_length, free_length_vertical_member, grid_height_m, Full_cross_member_num, Full_cross_member_length,
              #Full_cross_member_total_length, uncomplete_cross_member_a_length, Full_diagonal_member_num, Full_diagonal_member_length, Full_diagonal_member_total_length, area_length_CtoC, p, area_width_CtoC, uncomplete_cross_member_a_num):
    fig, ax = plt.subplots(figsize=(8.27, 11.69))
     # Create a new DXF document
    doc = ezdxf.new()
    msp = doc.modelspace()  # Model space for adding elements

    for i in range(num_main_members):
        y = free_length_vertical_member + i * grid_height_m
        ax.plot([0, area_length], [y, y], color='orange', label='Main Member Line' if i == 0 else "")
        msp.add_line((0, y*1000), (area_length*1000, y*1000    ))  # Save line to DXF

    # Plot uncomplete cross members top left
    slope = math.tan(math.radians(60))
    x_start = (area_length - area_length_CtoC)/2
    x_end = (area_length - area_length_CtoC)/2 + grid_size_m  # x_end = (p*math.tan(math.radians(30)) + grid_size_m*0.25) + (area_length -area_length_CtoC)/2
    y_end = area_width - (area_width -area_width_CtoC)/2
    y_start = slope*(x_start - x_end) + y_end
    first_iteration = True
    while y_start > (area_width - area_width_CtoC) / 2 and x_end < area_length:   
      # Plot the line with the current x_start and y_start values
      ax.plot([x_start, x_end],
                [y_start, y_end], color='green', label='Cross Member Line' if first_iteration else "")
      msp.add_line((x_start*1000, y_start*1000), (x_end*1000, y_end*1000))  # Save line to DXF

      y_start = slope*(x_start - x_end) + y_end

      # Update x_end and y_start for the next iteration
      x_end += grid_size_m  # Move right by grid size
      y_start = slope * (x_start - x_end) + y_end  # Recalculate y_start based on the new x_end
      first_iteration = False

    

    # Plot full cross members at a 60-degree angle
    if area_length >= area_width / (3 ** 0.5):
        slope = math.tan(math.radians(60))
        y_end = area_width - (area_width -area_width_CtoC)/2
        y_start = (area_width - area_width_CtoC)/2 #slope*(x_start - x_end) + y_end
        x_end = (area_length - area_length_CtoC)/2 + (uncomplete_cross_member_a_num + 1) * grid_size_m       #(p*math.tan(math.radians(30)) + grid_size_m*0.25) + (area_length -area_length_CtoC)/2 + (uncomplete_cross_member_a_num) * grid_size_m
        x_start = (y_start - y_end)/slope + x_end
        first_iteration = True
        while x_end < area_length - (area_length - area_length_CtoC)/2 :
          # Plot the line with the current x_start and y_start values
          ax.plot([x_start, x_end], [y_start, y_end], color='green') #, label='Full Cross Member' if first_iteration else "")
          msp.add_line((x_start*1000, y_start*1000), (x_end*1000, y_end*1000))  # Save line to DXF
          x_end = x_end + grid_size_m
          x_start = (y_start - y_end)/slope + x_end
          first_iteration = False
    else:
        # Case when area_width > area_length / (3 ** 0.5)
        slope = math.tan(math.radians(60))

        x_start = (area_length - area_length_CtoC) / 2
        y_start = area_width - ((area_width - area_width_CtoC) / 2 + grid_height_m * 2 * (uncomplete_cross_member_a_num + 1))
        x_end = area_length - (area_length - area_length_CtoC) / 2
        y_end = slope * (x_end - x_start) + y_start

        first_iteration = True

        while y_start > (area_width - area_width_CtoC) / 2:
           ax.plot([x_start, x_end], [y_start, y_end], color='green') #, label='Cross Member Line' if first_iteration else "")
           msp.add_line((x_start * 1000, y_start * 1000), (x_end * 1000, y_end * 1000))  # Save line to DXF

           # Move downward by 2 * grid_height_m
           y_start -= 2 * grid_height_m
           y_end = slope * (x_end - x_start) + y_start  # Recalculate y_end

           first_iteration = False

    # Plot uncomplete cross members Bottom right
    slope = math.tan(math.radians(60))
    x_end = area_length - (area_length - area_length_CtoC)/2
    x_start = (((area_width - area_width_CtoC)/2) - (area_width - (area_width -area_width_CtoC)/2))/slope + (area_length - area_length_CtoC)/2 + (uncomplete_cross_member_a_num + 1) * grid_size_m + (Full_cross_member_num )*grid_size_m   
    y_start = (area_width -area_width_CtoC)/2 #slope*(x_start - x_end) + y_end
    y_end = slope*(x_end - x_start) + y_start
    first_iteration = True  # Flag to track first iteration

    while x_start < area_length - (area_length - area_length_CtoC)/2 :
     #Plot the line with the current x_start and y_start values
      #if y_end > 0:
      ax.plot([x_start, x_end], [y_start, y_end], color='green')
      msp.add_line((x_start*1000, y_start*1000), (x_end*1000, y_end*1000))  # Save line to DXF
      x_start = x_start + grid_size_m
      y_end = slope*(x_end - x_start) + y_start
      first_iteration = False



    # Plot uncomplete diagonal members Bottom left
    slope = - math.tan(math.radians(60))
    x_start = (area_length - area_length_CtoC)/2
    x_end = (((area_width - area_width_CtoC)/2) - (area_width - (area_width -area_width_CtoC)/2))/slope + ((area_length - area_length_CtoC)/2 + grid_size_m - 2*p*math.tan(math.radians(30))) - grid_size_m
     #(area_length - area_length_CtoC)/2 + (grid_size_m-2*p/math.tan(math.radians(60)))+ abs(grid_size_m - 2*p/math.cos(math.radians(30)))*0.5
    y_end = (area_width - area_width_CtoC)/2
    y_start = slope*(x_start - x_end) + y_end
    first_iteration = True
    while (area_length - area_length_CtoC)/2 <= x_end :
      
        if x_end < (area_length + area_length_CtoC)/2:
            # Plot the line with the current x_start and y_start values
            ax.plot([x_start, x_end],
                [y_start, y_end], color='cyan', label='Diagonal Member Line' if first_iteration else "")
            msp.add_line((x_start*1000, y_start*1000), (x_end*1000, y_end*1000))  # Save line to DXF
            first_iteration = False

        # Update x_end and y_start for the next iteration
        x_end -= grid_size_m  # Move right by grid size
        y_start = slope * (x_start - x_end) + y_end  # Recalculate y_start based on the new x_end
        

     # Plot full diagonal Member members at a 60-degree angle
    if area_length >= area_width / (3 ** 0.5):
        slope = math.tan(math.radians(120))
        y_start = area_width - (area_width -area_width_CtoC)/2
        y_end = (area_width - area_width_CtoC)/2 #slope*(x_start - x_end) + y_end
        #x_end = (area_length - area_length_CtoC)/2 + (uncomplete_cross_member_a_num + 1) * grid_size_m       #(p*math.tan(math.radians(30)) + grid_size_m*0.25) + (area_length -area_length_CtoC)/2 + (uncomplete_cross_member_a_num) * grid_size_m
        x_start = (area_length - area_length_CtoC)/2 + grid_size_m - 2*p*math.tan(math.radians(30)) #(y_start - y_end)/slope + x_end
        x_end = (y_end - y_start)/slope + x_start
        first_iteration = True
 
        while x_end < area_length - (area_length - area_length_CtoC)/2 :
            
            # Plot the line with the current x_start and y_start values
            ax.plot([x_start, x_end], [y_start, y_end], color='cyan') #, label='Full diagonal Member' if first_iteration else "")
            msp.add_line((x_start*1000, y_start*1000), (x_end*1000, y_end*1000))  # Save line to DXF
            x_start = x_start + grid_size_m
            x_end = (y_end - y_start)/slope + x_start
            
    else:
        slope = math.tan(math.radians(120))
        x_start = (area_length - area_length_CtoC)/2
        x_end = (area_length + area_length_CtoC)/2
        y_start = (num_main_members-1) * grid_height_m + (area_width -area_width_CtoC)/2
        y_end = y_start + (x_end-x_start)*slope

        while y_end >= (area_width -area_width_CtoC)/2:
            # Plot the line with the current x_start and y_start values
            ax.plot([x_start, x_end], [y_start, y_end], color='cyan') #, label='Full diagonal Member' if first_iteration else "")
            msp.add_line((x_start*1000, y_start*1000), (x_end*1000, y_end*1000))  # Save line to DXF
            y_start = y_start - 2 * grid_height_m
            y_end = y_start + (x_end-x_start)*slope        
            

        # Plot uncomplete diagonal Member top right
    slope = -math.tan(math.radians(60))
    x_end = area_length - (area_length - area_length_CtoC)/2
    x_start = area_length - (area_length_CtoC-(int((area_length_CtoC-(grid_size_m-2*p/math.tan(math.radians(60))))/grid_size_m))*grid_size_m-(grid_size_m-2*p/math.tan(math.radians(60))) + (area_length - area_length_CtoC)/2)
     #((area_length - area_length_CtoC)/2) + ((uncomplete_cross_member_a_num + Full_cross_member_num) * grid_size_m + (grid_size_m - 2 * p * math.tan(math.radians(30)))) - grid_size_m
    y_start = area_width_CtoC + (area_width -area_width_CtoC)/2 #slope*(x_start - x_end) + y_end
    y_end = slope*(x_end - x_start) + y_start
    first_iteration = True  # Flag to track first iteration

    while y_end > (area_width -area_width_CtoC)/2 and x_start >= (area_length - area_length_CtoC)/2:
      #Plot the line with the current x_start and y_start values
      ax.plot([x_start, x_end], [y_start, y_end], color='cyan')
      msp.add_line((x_start*1000, y_start*1000), (x_end*1000, y_end*1000))  # Save line to DXF
      x_start = x_start - grid_size_m
      y_end = slope*(x_end - x_start) + y_start
      first_iteration = False


    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(0, area_length_CtoC + 1)
    ax.set_ylim(0, area_width + 1)
    # ax.set_title("Grid Pattern of - Opening", {cell_opening_size}, "mm ×", {cell_opening_size}, "mm ×", {cell_opening_size}, "mm")
    ax.set_title(f"Grid Pattern of - Opening {cell_opening_size}mm × {cell_opening_size}mm × {cell_opening_size}mm")
    ax.set_xlabel(f"Length of Area (m): {area_length} m")
    ax.set_ylabel(f"Width of Area (m): {area_width} m")
    # Hide axis lines but keep labels
    ax.set_xticks([])  # Hide x-axis ticks
    ax.set_yticks([])  # Hide y-axis ticks

    # Hide the spines (the outer box)
    for spine in ax.spines.values():
        spine.set_visible(False)
    # Remove axes and grid
    #ax.axis('off')         # Hides x and y axes
    plt.grid(False)        # Disables grid lines
    legend_width = 0.5 * area_width
    legend_height = 0.5 * area_length
    # Transform data coordinates to figure coordinates
    transform = ax.transData.transform
    inv_transform = ax.transAxes.inverted().transform

    # Convert (area_length, area_width) from data coordinates to figure coordinates
    legend_x, legend_y = inv_transform(transform((area_length, area_width)))
    # Add legend with bottom-left corner at (area_length, area_width)
    ax.legend(loc='lower left', bbox_to_anchor=(legend_x, legend_y), bbox_transform=ax.transAxes, fontsize=8, frameon=True)

    #ax.legend(loc='upper right')
    plt.grid()
    doc.saveas("grid_plot.dxf")
    plt.show()

    # Define a safer directory path, e.g., within the current working directory
    output_dir = "./output"

    # Ensure the directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Set the full paths for saving
    triangular_open_pdf_path = os.path.join(output_dir, "grid_plot.pdf")
    triangular_open_image_path = os.path.join(output_dir, "grid_plot.png")
    triangular_open_dxf_path = os.path.join(output_dir, "grid_plot.dxf")
    #triangular_open_cross_excel_path = os.path.join(output_dir, "grid_data.xlsx")

    # Save the figure as a PDF and PNG
    plt.savefig(triangular_open_pdf_path, format="pdf", bbox_inches="tight")  # Save as PDF
    plt.savefig(triangular_open_image_path, format="png", bbox_inches="tight")  # Save as PNG
    doc.saveas(triangular_open_dxf_path)  # Save DXF file

    # Save Excel to BytesIO (for Gradio download)
    # Save Excel to BytesIO (for Gradio download)
    excel_output = BytesIO()
    with pd.ExcelWriter(excel_output, engine="openpyxl") as writer:
        main_member_df.to_excel(writer, index=False, sheet_name=f"Main Member- {total_main_member} nos")
        cross_member_df.to_excel(writer, index=False, sheet_name=f"Cross Members- {total_cross_member} nos")
        diagonal_member_df.to_excel(writer, index=False, sheet_name=f"Diagonal Member- {total_diagonal_member} nos")
    excel_output.seek(0)

    main_member_df = main_member_df.fillna('')
    cross_member_df = cross_member_df.fillna('')
    diagonal_member_df = diagonal_member_df.fillna('')

    # Generate HTML tables for display
    html_output = f"""
    <div style='font-family: Arial, sans-serif;'>
        <h2>Main Members Details: {total_main_member} nos</h2>
        {main_member_df.to_html(index=False, classes='dataframe')}

        <h2>Cross Members Details: {total_cross_member} nos</h2>
        {cross_member_df.to_html(index=False, classes='dataframe')}

        <h2>Diagonal Members Details: : {total_diagonal_member} nos</h2>
        {diagonal_member_df.to_html(index=False, classes='dataframe')}
        
        <style>
            .dataframe {{
                width: 100%;
                border-collapse: collapse;
                margin-bottom: 20px;
                color: #000000;
                background-color: #ffffff;
            }}
            .dataframe th {{
                background-color: #ffa8B6;
                color: #373435;
                padding: 10px;
                text-align: center;
            }}
            .dataframe td {{
                padding: 8px;
                border-bottom: 1px solid #ddd;
                text-align: center;
            }}
            .dataframe tr:nth-child(even) {{
                background-color: #f2f2f2;
            }}

            /* Dark theme overrides */
            @media (prefers-color-scheme: dark) {{
                .dataframe {{
                    color: #ffffff;
                    background-color: #121212;
                }}
                .dataframe th {{
                    background-color: #444;
                    color: #ffa8b6;
                }}
                .dataframe td {{
                    border-bottom: 1px solid #444;
                }}
                .dataframe tr:nth-child(even) {{
                    background-color: #1e1e1e;
                }}
            }}
        </style>
    </div>
    """

    # Save to disk for Gradio
    excel_file_path = os.path.join(output_dir, "triangular_grid_output.xlsx")
    with open(excel_file_path, "wb") as f:
        f.write(excel_output.getvalue())
    
    # ✅ Final return
    return (
        triangular_open_grid_details,
        triangular_open_image_path,
        html_output,
        excel_file_path,
        triangular_open_pdf_path,
        triangular_open_dxf_path
    )

def draw_triangular_close_ended(area_length, area_width, cell_opening_size, baffle_width):
    # Convert grid dimensions from mm to meters
    cell_opening_size_m = cell_opening_size / 1000  # mm to meters
    grid_size_m = (cell_opening_size + baffle_width * math.tan(math.radians(60))) / 1000
    grid_height_m = grid_size_m * math.cos(math.radians(30))
    area_length_CtoC = area_length - (baffle_width * math.cos(math.radians(30))) / 1000
    area_width_CtoC = area_width - (baffle_width * math.cos(math.radians(60))) / 1000

    triangular_close_image_path = ""
    triangular_close_pdf_path = ""
    triangular_close_dxf_path = ""

    # Validate that grid dimensions do not exceed area dimensions
    if cell_opening_size_m >= area_width or cell_opening_size_m >= area_length:
        return "Grid dimensions must be less than or equal to area dimensions." , None, None

    # Validate that baffle width does not exceed grid dimensions
    if baffle_width >= cell_opening_size:
        return "Baffle width must be less than both grid length and width." , None, None

    # Calculate number of main members
    num_main_members = int(area_width / grid_height_m)   # Horizontal lines

    # Calculate initial free lengths
    free_length_vertical_member = (area_width - (num_main_members - 1) * grid_height_m)

    # Adjust number of cross members based on free length of main members
    free_length_vertical_member_old = free_length_vertical_member

    # if 0 <= free_length_vertical_member_old < 0.10 * grid_height_m:
    #     num_main_members -= 1

    # elif 0.10 * grid_height_m <= free_length_vertical_member_old < 0.20 * grid_height_m:
    #     if 0.60 * grid_height_m > 0.4:
    #         pass
    #     else:
    #         num_main_members -= 1
  
    # elif 0.20 * grid_height_m <= free_length_vertical_member_old < 0.60 * grid_height_m:
    #     pass

    # elif 0.60 * grid_height_m <= free_length_vertical_member_old < 0.70 * grid_height_m:
    #     if 0.60 * grid_height_m > 0.4:
    #         num_main_members += 1
  
    if 0.80 * grid_height_m <= free_length_vertical_member_old < 1.00 * grid_height_m:
         num_main_members += 1

    # Recalculate the free lengths with the adjusted number of main members
    free_length_vertical_member = (area_width - (num_main_members - 1) * grid_height_m) / 2

    # Calculate total lengths after adjustments
    total_main_length = num_main_members * area_length

    #Full_cross_member_length = (area_width * 1000 - baffle_width * math.cos(math.radians(60)))/(math.cos(math.radians(30)) * 1000)
    x_f = min(area_length_CtoC, area_width_CtoC / math.tan(math.radians(60)))
    
    # Calculate y_f using the equation y_f = tan(30°) * x_f
    y_f = math.tan(math.radians(60)) * x_f

    # Compute the distance between A(0,0) and B(x_f, y_f)
    Full_cross_member_length = math.sqrt(x_f**2 + y_f**2)

    p = (area_width_CtoC - (num_main_members - 1) * grid_size_m * math.cos(math.radians(30))) / 2

    i = 0
    standard_length = 0
    while True:
        next_length = grid_size_m * (i + 1) * 1000
        if next_length > 3000:
            break
        
        i += 1
        standard_length = grid_size_m * i * 1000

    # Total Length of uncomplete cross member_a
    # First intersection will be (grid_size_m)/2 far form the center width line on top main runner
    uncomplete_cross_member_a_length = 0  # Initialize total length
    i = 0
    j = 1
    # Initialize output string

    uncomplete_cross_member_a_num  = 0  # Initialize a uncomplete_cross_member_a_num for non-zero values
    total_cross_member = 0
    n = int(Full_cross_member_length/3)
    # Create DataFrames for all three tables
    cross_member_df = pd.DataFrame(columns=[
        "Cross Member line", "Total Length (in mm)"
        ] + [
            item for k in range(n+1) for item in [
                f"Member Name {k+1}",
                f"Length {k+1} (in mm)",
                f"First Punch {k+1} (in mm)"
            ]
        ])
    
    c_1 = 2*grid_size_m*( 1 + i) # c_1 = 2 * (p*math.tan(math.radians(30)) + grid_size_m*0.25 - (area_length -area_length_CtoC)/2 + grid_size_m * (i))  # Initialize the first value for c
    while c_1 < Full_cross_member_length:  # Loop until c exceeds the full length
        i += 1  # Increment i
        uncomplete_cross_member_a_length += c_1  # Add c to the total length
        if c_1 > 0:  # If c is non-zero, increment the uncomplete_cross_member_a_num
             #triangular_open_grid_details += f"- C_{j}: {c_1*1000:.2f} mm \n"
             uncomplete_cross_member_a_num += 1
             new_row = pd.DataFrame([[f"C{j}", f"{c_1 * 1000:.2f}"]], columns=["Cross Member line", "Total Length (in mm)"])
             cross_member_df = pd.concat([cross_member_df, new_row], ignore_index=True)
             total_length = c_1 * 1000
             solution_found = False
             tolerance = 0.01  # Acceptable rounding error

             if total_length < 3000:                    
                cross_member_df.at[cross_member_df.index[-1], f"Member Name {1}"] = f"C{j}_1"
                cross_member_df.at[cross_member_df.index[-1], f"Length {1} (in mm)"] = f"{c_1 * 1000:.2f}"
                cross_member_df.at[cross_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(p * 1000 - baffle_width + baffle_width/4) / math.cos(math.radians(30)):.2f}"
             else:
                A = grid_size_m * 500
                B = p * 1000 / math.cos(math.radians(30))
                C = 1.5 * grid_size_m * 1000

                y_max = int((3000 - A - B) // (grid_size_m * 1000))

                for y in reversed(range(y_max + 1)):
                    first_part_length = A + B + y * grid_size_m * 1000 

                    z_max = int((3000 - (C - B)) // (grid_size_m * 1000))
                    for z in reversed(range(z_max + 1)):
                        last_part_length = C - B + z * grid_size_m * 1000
                        remaining = total_length - first_part_length - last_part_length
                        n = remaining / standard_length

                        if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                            n = round(n)
                            total_check = first_part_length + n * standard_length + last_part_length
                            if abs(total_length - total_check) <= tolerance:
                                cross_member_df.at[cross_member_df.index[-1], f"Member Name {1}"] = f"C{j}_1"
                                cross_member_df.at[cross_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                                cross_member_df.at[cross_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(cell_opening_size/2):.2f}"
                                
                                # Add last part (comes after n standard parts + first)
                                cross_member_df.at[cross_member_df.index[-1], f"Member Name {n+2}"] = f"C{j}_{n+2}"
                                cross_member_df.at[cross_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                                cross_member_df.at[cross_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_opening_size/2):.2f}"      
                                for q in range(n):
                                    cross_member_df.at[cross_member_df.index[-1], f"Member Name {q+2}"] = f"C{j}_{q+2}"
                                    cross_member_df.at[cross_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                    cross_member_df.at[cross_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_opening_size/ 2):.2f}"                            
                                                                        
                                total_cross_member = total_cross_member + n + 2
                                solution_found = True                                    
                                break

                    if solution_found:
                        break

                if not solution_found:
                    print("No valid combination found.")

             j += 1
        c_1 = 2*grid_size_m*( 1 + i)  # Update c for the next iteration    

    Full_cross_member_num = max(int(((area_length - baffle_width * (math.cos(math.radians(30))) / 1000) - uncomplete_cross_member_a_num * grid_size_m) / grid_size_m),int((area_width_CtoC /(2 * grid_height_m))) - uncomplete_cross_member_a_num)
    #Full_cross_member_num = max(int(((area_length - baffle_width * (math.cos(math.radians(30))) / 1000) - uncomplete_cross_member_a_num * grid_size_m) / grid_size_m),int((area_width_CtoC - 2 * grid_height_m*uncomplete_cross_member_a_num ) / 2 * grid_height_m ))
    Full_cross_member_total_length = Full_cross_member_length * Full_cross_member_num
    Right_free_length = (area_length_CtoC - (grid_size_m - p / math.tan(math.radians(60)))) - (int((area_length_CtoC - (grid_size_m - p / math.tan(math.radians(60))))/grid_size_m))*grid_size_m
    
    if area_length >= area_width / (3 ** 0.5): 
        for i in range(Full_cross_member_num):
            new_row = pd.DataFrame(
                [[f"C{i + 1 + uncomplete_cross_member_a_num}", f"{Full_cross_member_length * 1000:.2f}"]],
                columns=["Cross Member line", "Total Length (in mm)"]
            )
            cross_member_df = pd.concat([cross_member_df, new_row], ignore_index=True)
            total_length = Full_cross_member_length * 1000
            solution_found = False
            tolerance = 0.01  # Acceptable rounding error  

            if total_length < 3000:                    
                cross_member_df.at[cross_member_df.index[-1], f"Member Name {1}"] = f"C{i + 1 + uncomplete_cross_member_a_num}_1"
                cross_member_df.at[cross_member_df.index[-1], f"Length {1} (in mm)"] = f"{total_length:.2f}"
                cross_member_df.at[cross_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(p * 1000 - baffle_width + baffle_width/4) / math.cos(math.radians(30)):.2f}"
            else:
                A = grid_size_m * 500
                B = p * 1000 / math.cos(math.radians(30))
                C = 1.5 * grid_size_m * 1000

                y_max = int((3000 - A - B) // (grid_size_m * 1000))
                
                for y in reversed(range(y_max + 1)):
                    first_part_length = A + B + y * grid_size_m * 1000                        

                    z_max = int((3000 - A - B) // (grid_size_m * 1000))
                    for z in reversed(range(z_max + 1)):
                        last_part_length = A + B + z * grid_size_m * 1000
                        remaining = total_length - first_part_length - last_part_length
                        n = remaining / standard_length
                        
                        if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                            n = round(n)
                            total_check = first_part_length + n * standard_length + last_part_length
                            if abs(total_length - total_check) <= tolerance:
                                cross_member_df.at[cross_member_df.index[-1], f"Member Name {1}"] = f"C{i + 1 + uncomplete_cross_member_a_num}_1"
                                cross_member_df.at[cross_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                                cross_member_df.at[cross_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(cell_opening_size/2):.2f}"
                                
                                # Add last part (comes after n standard parts + first)
                                cross_member_df.at[cross_member_df.index[-1], f"Member Name {n+2}"] = f"C{i + 1 + uncomplete_cross_member_a_num}_{n+2}"
                                cross_member_df.at[cross_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                                cross_member_df.at[cross_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_opening_size/2):.2f}"      
                                for q in range(n):
                                    cross_member_df.at[cross_member_df.index[-1], f"Member Name {q+2}"] = f"C{i + 1 + uncomplete_cross_member_a_num}_{q+2}"
                                    cross_member_df.at[cross_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                    cross_member_df.at[cross_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_opening_size/ 2):.2f}"                                      
                                total_cross_member = total_cross_member + n + 2
                                solution_found = True                                    
                                break

                    if solution_found:
                        break

                if not solution_found:
                    print("No valid combination found.")
    else:
        for i in range(Full_cross_member_num):
            new_row = pd.DataFrame([[f"C{i + 1 + uncomplete_cross_member_a_num}", f"{Full_cross_member_length * 1000:.2f}"]], columns=["Cross Member line", "Total Length (in mm)"])
          
        
            cross_member_df = pd.concat([cross_member_df, new_row], ignore_index=True)
            total_length = Full_cross_member_length * 1000
            solution_found = False
            tolerance = 0.01  # Acceptable rounding error  

            if total_length < 3000:                    
                cross_member_df.at[cross_member_df.index[-1], f"Member Name {1}"] = f"C{i + 1 + uncomplete_cross_member_a_num}_1"
                cross_member_df.at[cross_member_df.index[-1], f"Length {1} (in mm)"] = f"{total_length:.2f}"
                cross_member_df.at[cross_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(p * 1000 - baffle_width + baffle_width/4) / math.cos(math.radians(30)):.2f}"
            else:
                A = grid_size_m * 500
                B = p * 1000 / math.cos(math.radians(30))
                C = 1.5 * grid_size_m * 1000


                y_max = int((3000 - A - 2*Right_free_length*1000) // (grid_size_m * 1000))

                for y in reversed(range(y_max + 1)):
                    first_part_length = 2*Right_free_length*1000 + A + y * grid_size_m * 1000                        

                    z_max = int((3000 - (C - B)) // (grid_size_m * 1000))
                    for z in reversed(range(z_max + 1)):
                        last_part_length = C - B + z * grid_size_m * 1000
                        remaining = total_length - first_part_length - last_part_length
                        n = remaining / standard_length

                        if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                            n = round(n)
                            total_check = first_part_length + n * standard_length + last_part_length
                            if abs(total_length - total_check) <= tolerance:
                                cross_member_df.at[cross_member_df.index[-1], f"Member Name {1}"] = f"C{i + 1 + uncomplete_cross_member_a_num}_1"
                                cross_member_df.at[cross_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                                cross_member_df.at[cross_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(cell_opening_size/2):.2f}"
                                
                                # Add last part (comes after n standard parts + first)
                                cross_member_df.at[cross_member_df.index[-1], f"Member Name {n+2}"] = f"C{i + 1 + uncomplete_cross_member_a_num}_{n+2}"
                                cross_member_df.at[cross_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                                cross_member_df.at[cross_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_opening_size/2):.2f}"      
                                for q in range(n):
                                    cross_member_df.at[cross_member_df.index[-1], f"Member Name {q+2}"] = f"C{i + 1 + uncomplete_cross_member_a_num}_{q+2}"
                                    cross_member_df.at[cross_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                    cross_member_df.at[cross_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_opening_size/ 2):.2f}"                                                               
                                total_cross_member = total_cross_member + n + 2
                                solution_found = True                                    
                                break

                    if solution_found:
                        break

                if not solution_found:
                    print("No valid combination found.")

    # Total Length of uncomplete cross member_b Bottom right
    uncomplete_cross_member_b_length = 0  # Initialize total length
    i = 0
    j = 1
    # triangular_open_grid_details += f"""
    # <b>Uncomplete Cross Member (Bottom Right) Details: </b><br>
    # """
    slope =  math.tan(math.radians(60))
    uncomplete_cross_member_b_num = 0  # Initialize a uncomplete_cross_member_b_num for non-zero values
    c_2 = ((area_length_CtoC - (((((area_width - area_width_CtoC)/2) - (area_width - (area_width -area_width_CtoC)/2))/slope + (area_length - area_length_CtoC)/2 + (uncomplete_cross_member_a_num + 1) * grid_size_m + (Full_cross_member_num - 1)*grid_size_m) - (area_length - area_length_CtoC)/2)) - i * grid_size_m)*2
    # (area_length_CtoC - ((uncomplete_cross_member_a_num + Full_cross_member_num) * grid_size_m + 0.5 * abs(grid_size_m - 2 * p / math.cos(math.radians(30)))) - grid_size_m + i * grid_size_m) * 2  # Initialize the first value for c  # Initialize the first value for c
      # (area_length_CtoC - ((uncomplete_cross_member_a_num + Full_cross_member_num) * grid_size_m + (grid_size_m - 2 * p * math.tan(math.radians(30)))) - grid_size_m + i * grid_size_m) * 2                #

    while c_2 > (area_length - area_length_CtoC)/2:  # Loop until c exceeds the full length
        i += 1  # Increment i

        if 0 < c_2 < Full_cross_member_length:  # If c is non-zero, increment the uncomplete_cross_member_b_num
            #triangular_open_grid_details += f"- C_{j+Full_cross_member_num+uncomplete_cross_member_a_num}: {c_2*1000:.2f} mm \n"
            uncomplete_cross_member_b_length += c_2  # Add c to the total length
            uncomplete_cross_member_b_num += 1  # Increment the counter
            new_row = pd.DataFrame([[f"C{j + Full_cross_member_num + uncomplete_cross_member_a_num}", f"{c_2 * 1000:.2f}"]], columns=["Cross Member line", "Total Length (in mm)"])
            cross_member_df = pd.concat([cross_member_df, new_row], ignore_index=True)
            total_length = c_2 * 1000
            solution_found = False
            tolerance = 0.01  # Acceptable rounding error

            if total_length < 3000:                    
                cross_member_df.at[cross_member_df.index[-1], f"Member Name {1}"] = f"C{j + Full_cross_member_num + uncomplete_cross_member_a_num}_1"
                cross_member_df.at[cross_member_df.index[-1], f"Length {1} (in mm)"] = f"{c_2 * 1000:.2f}"
                cross_member_df.at[cross_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(p * 1000 - baffle_width + baffle_width/4) / math.cos(math.radians(30)):.2f}"
            else:
                A = grid_size_m * 500
                B = p * 1000 / math.cos(math.radians(30))
                C = 1.5 * grid_size_m * 1000    

                y_max = int((3000 - A - Right_free_length*2000) // (grid_size_m * 1000))

                for y in reversed(range(y_max + 1)):
                    first_part_length = A + Right_free_length*2000 + y * grid_size_m * 1000                        

                    z_max = int((3000 - (A + B)) // (grid_size_m * 1000))
                    for z in reversed(range(z_max + 1)):
                        last_part_length = A +  B + z * grid_size_m * 1000
                        remaining = total_length - first_part_length - last_part_length
                        n = remaining / standard_length

                        if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                            n = round(n)
                            total_check = first_part_length + n * standard_length + last_part_length
                            if abs(total_length - total_check) <= tolerance:
                                cross_member_df.at[cross_member_df.index[-1], f"Member Name {1}"] = f"C{j + Full_cross_member_num + uncomplete_cross_member_a_num}_1"
                                cross_member_df.at[cross_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                                cross_member_df.at[cross_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(cell_opening_size/2):.2f}"
                                
                                # Add last part (comes after n standard parts + first)
                                cross_member_df.at[cross_member_df.index[-1], f"Member Name {n+2}"] = f"C{j + Full_cross_member_num + uncomplete_cross_member_a_num}_{n+2}"
                                cross_member_df.at[cross_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                                cross_member_df.at[cross_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_opening_size/2):.2f}"      
                                for q in range(n):
                                    cross_member_df.at[cross_member_df.index[-1], f"Member Name {q+2}"] = f"C{j + Full_cross_member_num + uncomplete_cross_member_a_num}_{q+2}"
                                    cross_member_df.at[cross_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                    cross_member_df.at[cross_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_opening_size/ 2):.2f}"                            
                                                                        
                                total_cross_member = total_cross_member + n + 2
                                solution_found = True                                    
                                break

                    if solution_found:
                        break
                if not solution_found:
                    print("No valid combination found.")
            j += 1
        c_2 =  ((area_length_CtoC - (((((area_width - area_width_CtoC)/2) - (area_width - (area_width -area_width_CtoC)/2))/slope + (area_length - area_length_CtoC)/2 + (uncomplete_cross_member_a_num + 1) * grid_size_m + (Full_cross_member_num - 1)*grid_size_m) - (area_length - area_length_CtoC)/2)) - i * grid_size_m)*2  # Update c for the next iteration

    
    n = int(Full_cross_member_length/3)
    diagonal_member_df = pd.DataFrame(columns=[
        "Diagonal Member line", "Total Length (in mm)"
        ] + [
            item for k in range(n+1) for item in [
                f"Member Name {k+1}",
                f"Length {k+1} (in mm)",
                f"First Punch {k+1} (in mm)"
            ]
        ])
    
    Full_diagonal_member_length = Full_cross_member_length #(area_width * 1000 - baffle_width * math.cos(math.radians(60))) / (math.cos(math.radians(30)) * 1000)
    total_diagonal_member = 0

    # Total Length of uncomplete diagonal member_a Top Left
    uncomplete_diagonal_member_a_length = 0  # Initialize total length
    i = 0
    j = 1
    # triangular_open_grid_details += f"""
    # <b>Uncomplete Diagonal Member (Top Right) Details: </b><br>
    # """    
    uncomplete_diagonal_member_a_num = 0  # Initialize a uncomplete_diagonal_member_a_num for non-zero values
    c_3 =  2*(area_length_CtoC-(int((area_length_CtoC-(grid_size_m-2*p/math.tan(math.radians(60))))/grid_size_m))*grid_size_m-(grid_size_m-2*p/math.tan(math.radians(60))) + i*grid_size_m)  # Initialize the first value for c
    while c_3 < Full_diagonal_member_length:  # Loop until c exceeds the full length
        i += 1  # Increment i


        if c_3 > 0:  # If c is non-zero, increment the uncomplete_diagonal_member_a_num
            #triangular_open_grid_details += f"- D_{j}: {c_3*1000:.2f} mm \n"
            uncomplete_diagonal_member_a_num += 1
            uncomplete_diagonal_member_a_length += c_3  # Add c to the total length
            new_row = pd.DataFrame([[f"D{j}", f"{c_3 * 1000:.2f}"]], columns=["Diagonal Member line", "Total Length (in mm)"])
            diagonal_member_df = pd.concat([diagonal_member_df, new_row], ignore_index=True)
            total_length = c_3 * 1000
            solution_found = False
            tolerance = 0.01  # Acceptable rounding error

            if total_length < 3000:                    
                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {1}"] = f"D{j}_1"
                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {1} (in mm)"] = f"{c_3 * 1000:.2f}"
                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(p * 1000 - baffle_width + baffle_width/4) / math.cos(math.radians(30)):.2f}"
            else:
                A = grid_size_m * 500
                B = p * 1000 / math.cos(math.radians(30))
                C = 1.5 * grid_size_m * 1000

                y_max = int((3000 - A - B) // (grid_size_m * 1000))
                for y in reversed(range(y_max + 1)):
                    first_part_length = A + B + y * grid_size_m * 1000

                    z_max = int((3000 - A - Right_free_length*2000) // (grid_size_m * 1000))
                    for z in reversed(range(z_max + 1)):
                        last_part_length = A + Right_free_length*2000 + z * grid_size_m * 1000
                        
                        remaining = total_length - first_part_length - last_part_length
                        n = remaining / standard_length

                        if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                            n = round(n)
                            total_check = first_part_length + n * standard_length + last_part_length
                            if abs(total_length - total_check) <= tolerance:
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {1}"] = f"D{j}_1"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(cell_opening_size/2):.2f}"
                                
                                # Add last part (comes after n standard parts + first)
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {n+2}"] = f"D{j}_{n+2}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_opening_size/2):.2f}"      
                                for q in range(n):
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {q+2}"] = f"D{j}_{q+2}"
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_opening_size/ 2):.2f}"                            
                                total_diagonal_member = total_diagonal_member + n + 2                                        
                                solution_found = True                                    
                                break

                    if solution_found:
                        break

                if not solution_found:
                    print("No valid combination found.")
            j += 1
        c_3 = 2*(area_length_CtoC-(int((area_length_CtoC-(grid_size_m-2*p/math.tan(math.radians(60))))/grid_size_m))*grid_size_m-(grid_size_m-2*p/math.tan(math.radians(60))) + i*grid_size_m)  # Update c for the next iteration

    # Print the number of non-zero iterations and the total length after the loop
    # triangular_open_grid_details += f"""  
    # - Number of uncomplete diagonal members in top right side: {uncomplete_diagonal_member_a_num}
    # - Total Length of uncomplete diagonal member in top right side: {uncomplete_diagonal_member_a_length:.2f}
    # """
    Full_diagonal_member_num = max(int(((area_length - baffle_width * (math.cos(math.radians(30))) / 1000) - uncomplete_diagonal_member_a_num * grid_size_m -(area_length_CtoC-(int((area_length_CtoC-(grid_size_m-2*p/math.tan(math.radians(60))))/grid_size_m))*grid_size_m-(grid_size_m-2*p/math.tan(math.radians(60))))) / grid_size_m) + 1, 
                                   int((area_width_CtoC - (area_length_CtoC-(int((area_length_CtoC-(grid_size_m-2*p/math.tan(math.radians(60))))/grid_size_m))*grid_size_m-(grid_size_m-2*p/math.tan(math.radians(60))))) /(2 * grid_height_m)) - uncomplete_cross_member_a_num + 1) 
  

    
    #Full_diagonal_member_num = int(((area_length - baffle_width * (math.cos(math.radians(30))) / 1000) - uncomplete_diagonal_member_a_num * grid_size_m) / grid_size_m)
    #(int((area_length_CtoC-(grid_size_m-2*p/math.tan(math.radians(60))))/grid_size_m)) - uncomplete_diagonal_member_a_num    
    Full_diagonal_member_total_length = Full_diagonal_member_length * Full_diagonal_member_num

    if area_length >= area_width / (3 ** 0.5): 
        for i in range(Full_diagonal_member_num):
            new_row = pd.DataFrame(
                [[f"C{i + 1 + uncomplete_diagonal_member_a_num}", f"{Full_diagonal_member_length * 1000:.2f}"]],
                columns=["Diagonal Member line", "Total Length (in mm)"]
            )
            diagonal_member_df = pd.concat([diagonal_member_df, new_row], ignore_index=True)
            total_length = Full_diagonal_member_length * 1000
            solution_found = False
            tolerance = 0.01  # Acceptable rounding error  

            if total_length < 3000:                    
                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {1}"] = f"D{i + 1 + uncomplete_diagonal_member_a_num}_1"
                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {1} (in mm)"] = f"{total_length:.2f}"
                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(p * 1000 - baffle_width + baffle_width/4) / math.cos(math.radians(30)):.2f}"
            else:
                A = grid_size_m * 500
                B = p * 1000 / math.cos(math.radians(30))
                C = 1.5 * grid_size_m * 1000

                y_max = int((3000 - A - B) // (grid_size_m * 1000))
                
                for y in reversed(range(y_max + 1)):
                    first_part_length = A + B + y * grid_size_m * 1000                        

                    z_max = int((3000 - A - B) // (grid_size_m * 1000))
                    for z in reversed(range(z_max + 1)):
                        last_part_length = A + B + z * grid_size_m * 1000
                        remaining = total_length - first_part_length - last_part_length
                        n = remaining / standard_length
                        
                        if 0 <= last_part_length < 3000 and n >= 0 and abs(n - round(n)) <= 1e-6:
                            n = round(n)
                            total_check = first_part_length + n * standard_length + last_part_length
                            if abs(total_length - total_check) <= tolerance:
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {1}"] = f"D{i + 1 + uncomplete_diagonal_member_a_num}_1"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(cell_opening_size/2):.2f}"
                                
                                # Add last part (comes after n standard parts + first)
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {n+2}"] = f"D{i + 1 + uncomplete_diagonal_member_a_num}_{n+2}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_opening_size/2):.2f}"      
                                for q in range(n):
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {q+2}"] = f"D{i + 1 + uncomplete_diagonal_member_a_num}_{q+2}"
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_opening_size/ 2):.2f}"                                      
                                total_diagonal_member = total_diagonal_member + n + 2
                                solution_found = True                                    
                                break

                    if solution_found:
                        break

                if not solution_found:
                    print("No valid combination found.")
    else:
        for i in range(Full_diagonal_member_num):
            new_row = pd.DataFrame([[f"C{i + 1 + uncomplete_diagonal_member_a_num}", f"{Full_diagonal_member_length * 1000:.2f}"]], columns=["Diagonal Member line", "Total Length (in mm)"])
          
        
            diagonal_member_df = pd.concat([diagonal_member_df, new_row], ignore_index=True)
            total_length = Full_diagonal_member_length * 1000
            solution_found = False
            tolerance = 0.01  # Acceptable rounding error  

            if total_length < 3000:                    
                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {1}"] = f"D{i + 1 + uncomplete_diagonal_member_a_num}_1"
                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {1} (in mm)"] = f"{total_length:.2f}"
                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(p * 1000 - baffle_width + baffle_width/4) / math.cos(math.radians(30)):.2f}"
            else:
                A = grid_size_m * 500
                B = p * 1000 / math.cos(math.radians(30))
                C = 1.5 * grid_size_m * 1000


                y_max = int((3000 - (C - B)) // (grid_size_m * 1000))

                for y in reversed(range(y_max + 1)):
                    first_part_length = C - B + y * grid_size_m * 1000                        

                    z_max = int((3000 - A - 2*Right_free_length*1000) // (grid_size_m * 1000))
                    for z in reversed(range(z_max + 1)):
                        last_part_length = 2*Right_free_length*1000 + A + z * grid_size_m * 1000
                        remaining = total_length - first_part_length - last_part_length
                        n = remaining / standard_length

                        if 0 <= last_part_length < 3000 and n >= 0 and abs(n - round(n)) <= 1e-6:
                            n = round(n)
                            total_check = first_part_length + n * standard_length + last_part_length
                            if abs(total_length - total_check) <= tolerance:
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {1}"] = f"D{i + 1 + uncomplete_diagonal_member_a_num}_1"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(cell_opening_size/2):.2f}"
                                
                                # Add last part (comes after n standard parts + first)
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {n+2}"] = f"D{i + 1 + uncomplete_diagonal_member_a_num}_{n+2}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_opening_size/2):.2f}"      
                                for q in range(n):
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {q+2}"] = f"D{i + 1 + uncomplete_diagonal_member_a_num}_{q+2}"
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_opening_size/ 2):.2f}"                                                               
                                total_diagonal_member = total_diagonal_member + n + 2
                                solution_found = True                                    
                                break

                    if solution_found:
                        break

                if not solution_found:
                    print("No valid combination found.")





    # Total Length of uncomplete diagonal member_b Left Bottom
    uncomplete_diagonal_member_b_length = 0  # Initialize total length
    i = 0
    j = 1
    # triangular_open_grid_details += f"""
    # <b>Uncomplete Diagonal Member (Bottom left) Details: </b><br>
    # """ 
    slope = - math.tan(math.radians(60))
    uncomplete_diagonal_member_b_num = 0  # Initialize a uncomplete_diagonal_member_b_num for non-zero values

    c_4 = (((((area_width - area_width_CtoC)/2) - (area_width - (area_width -area_width_CtoC)/2))/slope + ((area_length - area_length_CtoC)/2 + grid_size_m - 2*p*math.tan(math.radians(30))) - grid_size_m) - (area_length - area_length_CtoC)/2 - i * grid_size_m)*2
     #((grid_size_m-2*p/math.tan(math.radians(60)))+abs(grid_size_m - 2*p/math.cos(math.radians(30)))*0.5 + i * grid_size_m)*2 # Initialize the first value for
     #(area_length_CtoC - ((uncomplete_cross_member_a_num + Full_cross_member_num) * grid_size_m + (grid_size_m - 2 * p * math.tan(math.radians(30)))) - grid_size_m + i * grid_size_m) * 2
    while c_4 > (area_length - area_length_CtoC)/2:  # Loop until c exceeds the full length
        i += 1  # Increment i

        if 0 < c_4 < Full_diagonal_member_length:  # If c_4 is non-zero, increment the uncomplete_diagonal_member_b_num
            #triangular_open_grid_details += f"- D_{j+Full_diagonal_member_num+uncomplete_diagonal_member_a_num}: {c_4*1000:.2f} mm\n"
            uncomplete_diagonal_member_b_length += c_4  # Add c to the total length
            uncomplete_diagonal_member_b_num += 1  # Increment the counter
            new_row = pd.DataFrame([[f"D{j + Full_diagonal_member_num + uncomplete_diagonal_member_a_num}", f"{c_4 * 1000:.2f}"]], columns=["Diagonal Member line", "Total Length (in mm)"])
            diagonal_member_df = pd.concat([diagonal_member_df, new_row], ignore_index=True)
            total_length = c_4 * 1000
            solution_found = False
            tolerance = 0.01  # Acceptable rounding error

            if total_length < 3000:                    
                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {1}"] = f"D{j + Full_diagonal_member_num + uncomplete_diagonal_member_a_num}_1"
                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {1} (in mm)"] = f"{c_4 * 1000:.2f}"
                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(p * 1000 - baffle_width + baffle_width/4) / math.cos(math.radians(30)):.2f}"
            else:
                A = grid_size_m * 500
                B = p * 1000 / math.cos(math.radians(30))
                C = 1.5 * grid_size_m * 1000    

                y_max = int((3000 - C + B) // (grid_size_m * 1000))

                for y in reversed(range(y_max + 1)):
                    first_part_length = C - B + y * grid_size_m * 1000                        

                    z_max = int((3000 - (A + B)) // (grid_size_m * 1000))
                    for z in reversed(range(z_max + 1)):
                        last_part_length = A +  B + z * grid_size_m * 1000
                        remaining = total_length - first_part_length - last_part_length
                        n = remaining / standard_length

                        if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                            n = round(n)
                            total_check = first_part_length + n * standard_length + last_part_length
                            if abs(total_length - total_check) <= tolerance:
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {1}"] = f"D{j + Full_diagonal_member_num + uncomplete_diagonal_member_a_num}_1"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(cell_opening_size/2):.2f}"
                                
                                # Add last part (comes after n standard parts + first)
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {n+2}"] = f"D{j + Full_diagonal_member_num + uncomplete_diagonal_member_a_num}_{n+2}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                                diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_opening_size/2):.2f}"      
                                for q in range(n):
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"Member Name {q+2}"] = f"D{j + Full_diagonal_member_num + uncomplete_diagonal_member_a_num}_{q+2}"
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                    diagonal_member_df.at[diagonal_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_opening_size/ 2):.2f}"                            
                                total_diagonal_member = total_diagonal_member + n + 2                                        
                                solution_found = True                                    
                                break

                    if solution_found:
                        break
                if not solution_found:
                    print("No valid combination found.")
            j += 1
        c_4 = (((((area_width - area_width_CtoC)/2) - (area_width - (area_width -area_width_CtoC)/2))/slope + ((area_length - area_length_CtoC)/2 + grid_size_m - 2*p*math.tan(math.radians(30))) - grid_size_m) - (area_length - area_length_CtoC)/2 - i * grid_size_m)*2  # Update c for the next iteration
        # (area_length_CtoC - ((uncomplete_cross_member_a_num + Full_cross_member_num) * grid_size_m + (grid_size_m - 2 * p * math.tan(math.radians(30)))) - grid_size_m + i * grid_size_m) * 2  #

    # Print the number of non-zero iterations and the total length after the loop
    # triangular_open_grid_details += f"""    
    # - Number of uncomplete diagonal members in Bottom left side: {uncomplete_diagonal_member_b_num}
    # - Total Length of uncomplete diagonal member in Bottom left side: {uncomplete_diagonal_member_b_length:.2f}
    # """

    Outer_member_total_length = 2*(area_length+area_width)

    Total_running_meter = total_main_length + Full_cross_member_total_length + uncomplete_cross_member_a_length + uncomplete_cross_member_b_length + Full_diagonal_member_total_length + uncomplete_diagonal_member_a_length + uncomplete_diagonal_member_b_length+Outer_member_total_length

    # ===== TABLE 3:Main Members Details ===== 
    n = int(area_length/3)
    main_member_df = pd.DataFrame(columns=[
        "Main Member line", "Total Length (in mm)"
        ] + [
            item for k in range(n+1) for item in [
                f"Member Name {k+1}",
                f"Length {k+1} (in mm)",
                f"First Punch {k+1} (in mm)"
            ]
        ])
    i = 1
    total_main_member = 0
    for i in range(num_main_members):
        if i % 2 == 0:
            new_row = pd.DataFrame(
                [[f"M{i + 1}", f"{area_length*1000:.2f}"]],
                columns=["Main Member line", "Total Length (in mm)"]
            )
            main_member_df = pd.concat([main_member_df, new_row], ignore_index=True)
            total_length = area_length * 1000
            solution_found = False
            tolerance = 0.01  # Acceptable rounding error  

            if total_length < 3000:                    
                main_member_df.at[main_member_df.index[-1], f"Member Name {1}"] = f"M{i + 1}_1"
                main_member_df.at[main_member_df.index[-1], f"Length {1} (in mm)"] = f"{total_length:.2f}"
                main_member_df.at[main_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(grid_size_m*1000 - (p * 1000 + 1.5 * baffle_width) * math.tan(math.radians(30)) + (area_length - area_length_CtoC)/2):.2f}"
            else:
                A = grid_size_m * 500
                B = p * 1000 / math.cos(math.radians(30))
                C = 1.5 * grid_size_m * 1000
                D = grid_size_m * 1000 - p * 1000 * math.tan(math.radians(30)) + (area_length - area_length_CtoC)*1000/2
                E = Right_free_length*1000 + (area_length - area_length_CtoC)*1000/2

                y_max = int((3000 - A - D) // (grid_size_m * 1000))
                
                for y in reversed(range(y_max + 1)):
                    first_part_length = A + D + y * grid_size_m * 1000 

                    z_max = int((3000 - A - E) // (grid_size_m * 1000))
                    for z in reversed(range(z_max + 1)):
                        last_part_length = A + E + z * grid_size_m * 1000
                        remaining = total_length - first_part_length - last_part_length
                        n = remaining / standard_length
                        
                        if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                            n = round(n)
                            total_check = first_part_length + n * standard_length + last_part_length
                            if abs(total_length - total_check) <= tolerance:
                                main_member_df.at[main_member_df.index[-1], f"Member Name {1}"] = f"M{i + 1}_1"
                                main_member_df.at[main_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                                main_member_df.at[main_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(cell_opening_size/2):.2f}"
                                
                                # Add last part (comes after n standard parts + first)
                                main_member_df.at[main_member_df.index[-1], f"Member Name {n+2}"] = f"M{i + 1}_{n+2}"
                                main_member_df.at[main_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                                main_member_df.at[main_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_opening_size/2):.2f}"      
                                for q in range(n):
                                    main_member_df.at[main_member_df.index[-1], f"Member Name {q+2}"] = f"M{i + 1}_{q+2}"
                                    main_member_df.at[main_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                    main_member_df.at[main_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_opening_size/ 2):.2f}"                                      
                                total_main_member = total_main_member + n + 2
                                solution_found = True                                    
                                break

                    if solution_found:
                        break

                if not solution_found:
                    print("No valid combination found.")
        else:
            new_row = pd.DataFrame(
                [[f"M{i + 1}", f"{area_length*1000:.2f}"]],
                columns=["Main Member line", "Total Length (in mm)"]
            )
            main_member_df = pd.concat([main_member_df, new_row], ignore_index=True)
            total_length = area_length * 1000
            solution_found = False
            tolerance = 0.01  # Acceptable rounding error  

            if total_length < 3000:                    
                main_member_df.at[main_member_df.index[-1], f"Member Name {1}"] = f"M{i + 1}_1"
                main_member_df.at[main_member_df.index[-1], f"Length {1} (in mm)"] = f"{total_length:.2f}"
                main_member_df.at[main_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(grid_size_m*1000 - (p * 1000 + 1.5 * baffle_width) * math.tan(math.radians(30)) + (area_length - area_length_CtoC)/2):.2f}"
            else:
                A = grid_size_m * 500
                B = p * 1000 / math.cos(math.radians(30))
                C = 1.5 * grid_size_m * 1000
                D = grid_size_m * 1000 - p * 1000 * math.tan(math.radians(30)) + (area_length - area_length_CtoC)*1000/2 - grid_size_m * 500
                E = Right_free_length*1000 + (area_length - area_length_CtoC)*1000/2 + grid_size_m * 500
                y_max = int((3000 - A - D) // (grid_size_m * 1000))
                
                for y in reversed(range(y_max + 1)):
                    first_part_length = A + D + y * grid_size_m * 1000 
                    z_max = int((3000 - A - E) // (grid_size_m * 1000))
                    for z in reversed(range(z_max + 1)):
                        last_part_length = A + E + z * grid_size_m * 1000

                        remaining = total_length - first_part_length - last_part_length
                        n = remaining / standard_length
                        
                        if 0 <= last_part_length < 3000 and n > -1 and abs(n - round(n)) <= 1e-6:
                            n = round(n)
                            total_check = first_part_length + n * standard_length + last_part_length
                            if abs(total_length - total_check) <= tolerance:
                                main_member_df.at[main_member_df.index[-1], f"Member Name {1}"] = f"M{i + 1}_1"
                                main_member_df.at[main_member_df.index[-1], f"Length {1} (in mm)"] = f"{first_part_length:.2f}"
                                main_member_df.at[main_member_df.index[-1], f"First Punch {1} (in mm)"] = f"{(cell_opening_size/2):.2f}"
                                
                                # Add last part (comes after n standard parts + first)
                                main_member_df.at[main_member_df.index[-1], f"Member Name {n+2}"] = f"M{i + 1}_{n+2}"
                                main_member_df.at[main_member_df.index[-1], f"Length {n+2} (in mm)"] = f"{last_part_length:.2f}"
                                main_member_df.at[main_member_df.index[-1], f"First Punch {n+2} (in mm)"] = f"{(cell_opening_size/2):.2f}"      
                                for q in range(n):
                                    main_member_df.at[main_member_df.index[-1], f"Member Name {q+2}"] = f"M{i + 1}_{q+2}"
                                    main_member_df.at[main_member_df.index[-1], f"Length {q+2} (in mm)"] = f"{standard_length:.2f}"
                                    main_member_df.at[main_member_df.index[-1], f"First Punch {q+2} (in mm)"] = f"{(cell_opening_size/ 2):.2f}"                                      
                                total_main_member = total_main_member + n + 2
                                solution_found = True                                    
                                break

                    if solution_found:
                        break

                if not solution_found:
                    print("No valid combination found.")


    remaining_outer_length_1 = 0
    remaining_outer_width_1 = 0

    full_outer_length_num = int((area_length-0.6)/3)
    remaining_outer_length = area_length*1000 - full_outer_length_num *3000 - 600
    if remaining_outer_length < 600 and remaining_outer_length > 0:
        full_outer_length_num -= 1
        remaining_outer_length_1 = area_length*1000 - full_outer_length_num *3000 - 600
    #print("remaining_outer_length =", remaining_outer_length)
    full_outer_width_num = int((area_width-0.6)/3)
    remaining_outer_width = area_width*1000 - full_outer_width_num *3000 - 600
    if remaining_outer_width < 600 and remaining_outer_width > 0:
        full_outer_width_num -= 1
        remaining_outer_width_1 = area_width*1000 - full_outer_width_num *3000 - 600

    #print("remaining_outer_width =", remaining_outer_width)

    total_outer_member = (full_outer_length_num + full_outer_width_num) * 2

    # Create outer_member_df
    outer_member_df = pd.DataFrame(columns=[
        "Product", "Quantity (NoS)"
    ])

    # Add full-length outer C Channel
    outer_member_df = pd.concat([
        outer_member_df,
        pd.DataFrame([["C Channel 25 × 103 × 3000 mm", f"{total_outer_member}"]], columns=["Product", "Quantity (NoS)"])
    ], ignore_index=True)

    # Check and add corner channels if necessary
    if 0 < remaining_outer_length < 600:
        outer_member_df = pd.concat([
            outer_member_df,
            pd.DataFrame([
                [f"C Channel 25 × 103  × {remaining_outer_length_1 / 2:.2f} mm", "4"]
            ], columns=["Product", "Quantity (NoS)"])
        ], ignore_index=True)

    elif remaining_outer_length >= 600:
        outer_member_df = pd.concat([
            outer_member_df,
            pd.DataFrame([
                [f"C Channel 25 × 103  × {remaining_outer_length:.2f} mm", "2"]
            ], columns=["Product", "Quantity (NoS)"])
        ], ignore_index=True)

    if 0 < remaining_outer_width < 600:
        outer_member_df = pd.concat([
            outer_member_df,
            pd.DataFrame([
                [f"C Channel 25 × 103 × 300 × {remaining_outer_width_1 / 2:.2f} mm", "4"]
            ], columns=["Product", "Quantity (NoS)"])
        ], ignore_index=True)

    elif remaining_outer_width >= 600:
        outer_member_df = pd.concat([
            outer_member_df,
            pd.DataFrame([
                [f"C Channel 25 × 103 × 300 × {remaining_outer_width:.2f} mm", "2"]
            ], columns=["Product", "Quantity (NoS)"])
        ], ignore_index=True)

    # If you always want to include the standard 300 mm corner channel regardless:
    outer_member_df = pd.concat([
        outer_member_df,
        pd.DataFrame([["Corner C Channel 25 × 103 × 300 × 300 mm", "4"]], columns=["Product", "Quantity (NoS)"])
    ], ignore_index=True)

    # Print results
    triangular_close_grid_details = "**Grid Details:**\n" 
    triangular_close_grid_details += f"""
    - triangular Grid Size Opening (in mm): {cell_opening_size} mm × {cell_opening_size} mm × {cell_opening_size} mm
    - Total required Running meter: {Total_running_meter:.2f} meters
    - Running Meter per Square Meter: {Total_running_meter/(area_length*area_width):.2f} running meters per square meter 
    - Main member Standard Punching Dimensions: {(cell_opening_size - 0.4):.2f} mm
    - Diaginal member Standard Punching Dimensions: {(cell_opening_size - 0.4):.2f} mm
    - Cross member Standard Punching Dimensions: {(grid_size_m*1000 - 29.5):.2f} mm
    """

    fig, ax = plt.subplots(figsize=(8.27, 11.69))
     # Create a new DXF document
    doc = ezdxf.new()
    msp = doc.modelspace()  # Model space for adding elements

    for i in range(num_main_members):
        y = free_length_vertical_member + i * grid_height_m
        ax.plot([0, area_length], [y, y], color='orange', label='Main Member' if i == 0 else "")
        msp.add_line((0, y*1000), (area_length*1000, y*1000    ))  # Save line to DXF

    # Plot uncomplete cross members top left
    slope = math.tan(math.radians(60))
    x_start = (area_length - area_length_CtoC)/2
    x_end = (area_length - area_length_CtoC)/2 + grid_size_m  # x_end = (p*math.tan(math.radians(30)) + grid_size_m*0.25) + (area_length -area_length_CtoC)/2
    y_end = area_width - (area_width -area_width_CtoC)/2
    y_start = slope*(x_start - x_end) + y_end
    first_iteration = True
    while y_start > (area_width - area_width_CtoC) / 2 and x_end < area_length:
      # Plot the line with the current x_start and y_start values
      ax.plot([x_start, x_end],
                [y_start, y_end], color='green', label='Cross Member Line' if first_iteration else "")
      msp.add_line((x_start*1000, y_start*1000), (x_end*1000, y_end*1000))  # Save line to DXF

      y_start = slope*(x_start - x_end) + y_end

      # Update x_end and y_start for the next iteration
      x_end += grid_size_m  # Move right by grid size
      y_start = slope * (x_start - x_end) + y_end  # Recalculate y_start based on the new x_end
      first_iteration = False

    # Plot full cross members at a 60-degree angle
    if area_length >= area_width/ (3 ** 0.5):
        slope = math.tan(math.radians(60))
        y_end = area_width - (area_width -area_width_CtoC)/2
        y_start = (area_width - area_width_CtoC)/2 #slope*(x_start - x_end) + y_end
        x_end = (area_length - area_length_CtoC)/2 + (uncomplete_cross_member_a_num + 1) * grid_size_m       #(p*math.tan(math.radians(30)) + grid_size_m*0.25) + (area_length -area_length_CtoC)/2 + (uncomplete_cross_member_a_num) * grid_size_m
        x_start = (y_start - y_end)/slope + x_end
        first_iteration = True
        while x_end < area_length - (area_length - area_length_CtoC)/2 :
          # Plot the line with the current x_start and y_start values
          ax.plot([x_start, x_end], [y_start, y_end], color='green') #, label='Full Cross Member' if first_iteration else "")
          msp.add_line((x_start*1000, y_start*1000), (x_end*1000, y_end*1000))  # Save line to DXF
          x_end = x_end + grid_size_m
          x_start = (y_start - y_end)/slope + x_end
          first_iteration = False
    else:
        # Case when area_width > area_length
        slope = math.tan(math.radians(60))

        x_start = (area_length - area_length_CtoC) / 2
        y_start = area_width - ((area_width - area_width_CtoC) / 2 + grid_height_m * 2 * (uncomplete_cross_member_a_num + 1))
        x_end = area_length - (area_length - area_length_CtoC) / 2
        y_end = slope * (x_end - x_start) + y_start

        first_iteration = True

        while y_start > (area_width - area_width_CtoC) / 2:
           ax.plot([x_start, x_end], [y_start, y_end], color='green') #, label='Cross Member Line' if first_iteration else "")
           msp.add_line((x_start * 1000, y_start * 1000), (x_end * 1000, y_end * 1000))  # Save line to DXF

           # Move downward by 2 * grid_height_m
           y_start -= 2 * grid_height_m
           y_end = slope * (x_end - x_start) + y_start  # Recalculate y_end

           first_iteration = False

    # Plot uncomplete cross members Bottom right
    slope = math.tan(math.radians(60))
    x_end = area_length - (area_length - area_length_CtoC)/2
    x_start = (((area_width - area_width_CtoC)/2) - (area_width - (area_width -area_width_CtoC)/2))/slope + (area_length - area_length_CtoC)/2 + (uncomplete_cross_member_a_num + 1) * grid_size_m + (Full_cross_member_num )*grid_size_m   
    y_start = (area_width -area_width_CtoC)/2 #slope*(x_start - x_end) + y_end
    y_end = slope*(x_end - x_start) + y_start
    first_iteration = True  # Flag to track first iteration

    while x_start < area_length - (area_length - area_length_CtoC)/2 :
      # Plot the line with the current x_start and y_start values
      #if y_end > 0:
      ax.plot([x_start, x_end], [y_start, y_end], color='green')
      msp.add_line((x_start*1000, y_start*1000), (x_end*1000, y_end*1000))  # Save line to DXF
      x_start = x_start + grid_size_m
      y_end = slope*(x_end - x_start) + y_start
      first_iteration = False



    # Plot uncomplete diagonal members Bottom left
    slope = - math.tan(math.radians(60))
    x_start = (area_length - area_length_CtoC)/2
    x_end = (((area_width - area_width_CtoC)/2) - (area_width - (area_width -area_width_CtoC)/2))/slope + ((area_length - area_length_CtoC)/2 + grid_size_m - 2*p*math.tan(math.radians(30))) - grid_size_m
     #(area_length - area_length_CtoC)/2 + (grid_size_m-2*p/math.tan(math.radians(60)))+ abs(grid_size_m - 2*p/math.cos(math.radians(30)))*0.5
    y_end = (area_width -area_width_CtoC)/2
    y_start = slope*(x_start - x_end) + y_end
    first_iteration = True
    while (area_length - area_length_CtoC)/2 <= x_end :
        if x_end < (area_length + area_length_CtoC)/2:
            # Plot the line with the current x_start and y_start values
            ax.plot([x_start, x_end],
                [y_start, y_end], color='cyan', label='Diagonal Member Line' if first_iteration else "")
            msp.add_line((x_start*1000, y_start*1000), (x_end*1000, y_end*1000))  # Save line to DXF
            first_iteration = False

        # Update x_end and y_start for the next iteration
        x_end -= grid_size_m  # Move right by grid size
        y_start = slope * (x_start - x_end) + y_end  # Recalculate y_start based on the new x_end
        

     # Plot full diagonal Member members at a 60-degree angle
    if area_length >= area_width/ (3 ** 0.5):
        slope = math.tan(math.radians(120))
        y_start = area_width - (area_width -area_width_CtoC)/2
        y_end = (area_width - area_width_CtoC)/2 #slope*(x_start - x_end) + y_end
        #x_end = (area_length - area_length_CtoC)/2 + (uncomplete_cross_member_a_num + 1) * grid_size_m       #(p*math.tan(math.radians(30)) + grid_size_m*0.25) + (area_length -area_length_CtoC)/2 + (uncomplete_cross_member_a_num) * grid_size_m
        x_start = (area_length - area_length_CtoC)/2 + grid_size_m - 2*p*math.tan(math.radians(30)) #(y_start - y_end)/slope + x_end
        x_end = (y_end - y_start)/slope + x_start
        first_iteration = True
 
        while x_end < area_length - (area_length - area_length_CtoC)/2 :
            
            # Plot the line with the current x_start and y_start values
            ax.plot([x_start, x_end], [y_start, y_end], color='cyan') #, label='Full diagonal Member' if first_iteration else "")
            msp.add_line((x_start*1000, y_start*1000), (x_end*1000, y_end*1000))  # Save line to DXF
            x_start = x_start + grid_size_m
            x_end = (y_end - y_start)/slope + x_start
            
    else:
        slope = math.tan(math.radians(120))
        x_start = (area_length - area_length_CtoC)/2
        x_end = (area_length + area_length_CtoC)/2
        y_start = (num_main_members-1) * grid_height_m + (area_width -area_width_CtoC)/2
        y_end = y_start + (x_end-x_start)*slope

        while y_end >= (area_width -area_width_CtoC)/2:
            # Plot the line with the current x_start and y_start values
            ax.plot([x_start, x_end], [y_start, y_end], color='cyan') #, label='Full diagonal Member' if first_iteration else "")
            msp.add_line((x_start*1000, y_start*1000), (x_end*1000, y_end*1000))  # Save line to DXF
            y_start = y_start - 2 * grid_height_m
            y_end = y_start + (x_end-x_start)*slope      


        # Plot uncomplete diagonal Member top right
    slope = -math.tan(math.radians(60))
    x_end = area_length - (area_length - area_length_CtoC)/2
    x_start = area_length - (area_length_CtoC-(int((area_length_CtoC-(grid_size_m-2*p/math.tan(math.radians(60))))/grid_size_m))*grid_size_m-(grid_size_m-2*p/math.tan(math.radians(60))) + (area_length - area_length_CtoC)/2)
     #((area_length - area_length_CtoC)/2) + ((uncomplete_cross_member_a_num + Full_cross_member_num) * grid_size_m + (grid_size_m - 2 * p * math.tan(math.radians(30)))) - grid_size_m
    y_start = area_width_CtoC + (area_width -area_width_CtoC)/2 #slope*(x_start - x_end) + y_end
    y_end = slope*(x_end - x_start) + y_start
    first_iteration = True  # Flag to track first iteration

    while y_end > (area_width - area_width_CtoC)/2 and x_start >= (area_length - area_length_CtoC)/2:
      #Plot the line with the current x_start and y_start values
      ax.plot([x_start, x_end], [y_start, y_end], color='cyan')
      msp.add_line((x_start*1000, y_start*1000), (x_end*1000, y_end*1000))  # Save line to DXF
      x_start = x_start - grid_size_m
      y_end = slope*(x_end - x_start) + y_start
      first_iteration = False

    # Plotting the additional lines
    ax.plot([0, area_length], [0, 0], color='blue', label='Outer Member')  # y = 0 line
    ax.plot([0, area_length], [area_width, area_width], color='blue') #, label='y=area_width' if num_main_members == 0 and num_cross_members == 0 else "")  # y = area_width line
    ax.plot([0, 0], [0, area_width], color='blue') #, label='x=0' if num_main_members == 0 and num_cross_members == 0 else "")  # x = 0 line
    ax.plot([area_length, area_length], [0, area_width], color='blue') #, label='x=area_length' if num_main_members == 0 and num_cross_members == 0 else "")  # x = area_length line

    # Create the outer members (DXF data)
    msp.add_line((0, 0), (area_length*1000, 0))  # x=0 line
    msp.add_line((0, area_width*1000), (area_length*1000, area_width*1000))  # x=area_length line
    msp.add_line((0, 0), (0, area_width*1000))  # y=0 line
    msp.add_line((area_length*1000, 0), (area_length*1000, area_width*1000))  # y=area_width line


    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(0, area_length_CtoC + 1)
    ax.set_ylim(0, area_width + 1)
    ax.set_title(f"Grid Pattern of - Opening {cell_opening_size}mm × {cell_opening_size}mm × {cell_opening_size}mm")
    ax.set_xlabel(f"Length of Area (m): {area_length} m")
    ax.set_ylabel(f"Width of Area (m): {area_width} m")
   # Hide axis lines but keep labels
    ax.set_xticks([])  # Hide x-axis ticks
    ax.set_yticks([])  # Hide y-axis ticks

    # Hide the spines (the outer box)
    for spine in ax.spines.values():
        spine.set_visible(False)
    # Remove axes and grid
    #ax.axis('off')         # Hides x and y axes
    plt.grid(False)        # Disables grid lines
    legend_width = 0.5 * area_width
    legend_height = 0.5 * area_length
    # Transform data coordinates to figure coordinates
    transform = ax.transData.transform
    inv_transform = ax.transAxes.inverted().transform

    # Convert (area_length, area_width) from data coordinates to figure coordinates
    legend_x, legend_y = inv_transform(transform((area_length, area_width)))
    # Add legend with bottom-left corner at (area_length, area_width)
    ax.legend(loc='lower left', bbox_to_anchor=(legend_x, legend_y), bbox_transform=ax.transAxes, fontsize=8, frameon=True)

    #ax.legend(loc='upper right')
    plt.grid()
    doc.saveas("grid_plot.dxf")
    plt.show()

    # Define a safer directory path, e.g., within the current working directory
    output_dir = "./output"

    # Ensure the directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Set the full paths for saving
    triangular_close_pdf_path = os.path.join(output_dir, "grid_plot.pdf")
    triangular_close_image_path = os.path.join(output_dir, "grid_plot.png")
    triangular_close_dxf_path = os.path.join(output_dir, "grid_plot.dxf")
    #triangular_close_cross_excel_path = os.path.join(output_dir, "grid_data.xlsx")

    # Save the figure as a PDF and PNG
    plt.savefig(triangular_close_pdf_path, format="pdf", bbox_inches="tight")  # Save as PDF
    plt.savefig(triangular_close_image_path, format="png", bbox_inches="tight")  # Save as PNG
    doc.saveas(triangular_close_dxf_path)  # Save DXF file

    # Save Excel to BytesIO (for Gradio download)
    # Save Excel to BytesIO (for Gradio download)
    excel_output = BytesIO()
    with pd.ExcelWriter(excel_output, engine="openpyxl") as writer:
        main_member_df.to_excel(writer, index=False, sheet_name=f"Main Memberr- {total_main_member} nos")
        cross_member_df.to_excel(writer, index=False, sheet_name=f"Cross Members- {total_cross_member} nos")
        diagonal_member_df.to_excel(writer, index=False, sheet_name=f"Diagonal Member- {total_diagonal_member} nos")
        outer_member_df.to_excel(writer, index=False, sheet_name="Outer Member")

    excel_output.seek(0)

    main_member_df = main_member_df.fillna('')
    cross_member_df = cross_member_df.fillna('')
    diagonal_member_df = diagonal_member_df.fillna('')
    outer_member_df = outer_member_df.fillna('')

    # Generate HTML tables for display
    html_output = f"""
    <div style='font-family: Arial, sans-serif;'>
        <h2>Main Members Details: {total_main_member} nos</h2>
        {main_member_df.to_html(index=False, classes='dataframe')}
    
        <h2>Cross Members Details: {total_cross_member} nos</h2>
        {cross_member_df.to_html(index=False, classes='dataframe')}
        
        <h2>Diagonal Members Details: {total_diagonal_member} nos</h2>
        {diagonal_member_df.to_html(index=False, classes='dataframe')}
        
        <h2>Outer Members Details</h2>
        {outer_member_df.to_html(index=False, classes='dataframe')}

        
        <style>
            .dataframe {{
                width: 100%;
                border-collapse: collapse;
                margin-bottom: 20px;
                color: #000000;
                background-color: #ffffff;
            }}
            .dataframe th {{
                background-color: #ffa8B6;
                color: #373435;
                padding: 10px;
                text-align: center;
            }}
            .dataframe td {{
                padding: 8px;
                border-bottom: 1px solid #ddd;
                text-align: center;
            }}
            .dataframe tr:nth-child(even) {{
                background-color: #f2f2f2;
            }}

            /* Dark theme overrides */
            @media (prefers-color-scheme: dark) {{
                .dataframe {{
                    color: #ffffff;
                    background-color: #121212;
                }}
                .dataframe th {{
                    background-color: #444;
                    color: #ffa8b6;
                }}
                .dataframe td {{
                    border-bottom: 1px solid #444;
                }}
                .dataframe tr:nth-child(even) {{
                    background-color: #1e1e1e;
                }}
            }}
        </style>
    </div>
    """

    # Save to disk for Gradio
    excel_file_path = os.path.join(output_dir, "triangular_grid_output.xlsx")
    with open(excel_file_path, "wb") as f:
        f.write(excel_output.getvalue())
    
    # ✅ Final return
    return (
        triangular_close_grid_details,
        triangular_close_image_path,
        html_output,
        excel_file_path,
        triangular_close_pdf_path,
        triangular_close_dxf_path
    )

import gradio as gr
import base64

# Function to dynamically update visible input fields based on the selected grid type
def change_inputs(grid_type):
    visibility_settings = {
        "Quadrilateral Open-Ended Grid": [True, True, True, True, False, True],
        "Quadrilateral Close-Ended Grid": [True, True, True, True, False, True],
        "Triangular Open-Ended Grid": [True, True, False, False, True, True],
        "Triangular Close-Ended Grid": [True, True, False, False, True, True]
    }
    return [gr.update(visible=visible) for visible in visibility_settings.get(grid_type, [False] * 6)]

# Function to handle drawing the grid
def draw_grid(grid_type, area_length, area_width, cell_length, cell_width, cell_opening_size, baffle_width):
    try:
        # Convert inputs to numbers, handling empty strings
        area_length = float(area_length) if area_length.strip() else 0
        area_width = float(area_width) if area_width.strip() else 0
        cell_length = float(cell_length) if cell_length.strip() else 0
        cell_width = float(cell_width) if cell_width.strip() else 0
        cell_opening_size = float(cell_opening_size) if cell_opening_size.strip() else 0
        #cell_ctoc_size = float(cell_ctoc_size) if cell_ctoc_size.strip() else 0
        baffle_width = float(baffle_width) if baffle_width.strip() else 0
    except ValueError:
        return "Error: Please enter valid numerical values.", None, None

    if grid_type == "Quadrilateral Open-Ended Grid":
        return draw_quadrilaterals_open_ended(area_length, area_width, cell_length, cell_width, baffle_width)
    elif grid_type == "Quadrilateral Close-Ended Grid":
        return draw_quadrilateral_close_ended(area_length, area_width, cell_length, cell_width, baffle_width)
    elif grid_type == "Triangular Open-Ended Grid":
        return draw_triangular_open_ended(area_length, area_width, cell_opening_size, baffle_width)
    elif grid_type == "Triangular Close-Ended Grid":
        return draw_triangular_close_ended(area_length, area_width, cell_opening_size, baffle_width)

# Function to clear inputs and outputs
def clear_inputs():
    return "", "", "", "", "", "", gr.update(visible=False), gr.update(visible=False), gr.update(visible=False), gr.update(visible=False), gr.update(visible=False)


html_code = """
<div style="display: grid; grid-template-columns: auto 1fr; align-items: center; gap: 16px; max-width: 10000px; margin: 0 auto;">
    <div style="grid-row: span 2; display: flex; align-items: center;">
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 595.3 841.9" style="height: 3.5em;">
            <style type="text/css">
                .st0 { fill: #71278F; }

                .st1 { fill: #373435; }
                @media (prefers-color-scheme: dark) {
                    .st1 { fill: #FFFFFF; }
                }
            </style>
            <g id="Design">
                <path class="st0" d="M0,222.2c0-5.8,5.8-11.7,11.7-11.7c40.9-5.8,70.2-17.5,99.4-35.1s46.8-40.9,64.3-70.2
                    c11.7-29.2,23.4-64.3,23.4-93.5c0-5.8,0-11.7,5.8-11.7h152c5.8,0,11.7,5.8,11.7,11.7v17.5c0,70.2-17.5,128.6-58.5,181.2
                    s-87.7,87.7-146.2,116.9c-58.5,23.4-105.2,35.1-152,35.1c-5.8,0-11.7,0-11.7-11.7V222.2z"/>
                <path class="st1" d="M169.5,374.2c40.9-17.5,81.9-40.9,116.9-64.3c35.1-29.2,64.3-64.3,87.7-105.2h146.2c5.8,0,11.7,5.8,11.7,11.7
                    v134.5c0,5.8-5.8,11.7-11.7,11.7H356.6v257.2c0,23.4,11.7,40.9,23.4,52.6c11.7,11.7,29.2,17.5,52.6,17.5h81.9
                    c5.8,0,11.7,5.8,11.7,11.7v128.6c0,5.8-5.8,11.7-11.7,11.7h-99.4c-76,0-134.5-17.5-175.4-58.5c-46.8-40.9-70.2-93.5-70.2-163.7
                    V374.2z"/>
            </g>
        </svg>
    </div>

    <h1 style="margin: 0 0 2px 0; grid-column: 2; font-size: 1.5em; line-height: 0.2;">BaffleWorks: Grid Layout & Cut Plan</h1>
    <p style="margin: 0; grid-column: 2; font-size: 0.95em; line-height: 0.01;">
        Enter the area dimensions, grid sizes, and baffle width to generate a grid layout with calculations.
        <a href="https://admintechno-my.sharepoint.com/:b:/g/personal/priyanshu_technoceiling_com/EcQlmnJC_7BAjxuFv7XFHwIBFfH9evlswYtM9By61iZstA" target="_blank">Link</a>.
    </p>
</div>
"""


# Gradio UI
with gr.Blocks() as interface:
    gr.HTML(html_code)
    grid_type = gr.Radio(
        ["Quadrilateral Open-Ended Grid", "Quadrilateral Close-Ended Grid", "Triangular Open-Ended Grid", "Triangular Close-Ended Grid"],
        label="Grid Type"
    )
    # Initialize the inputs, setting them visible=False initially
    area_length = gr.Textbox(label="Length of Area (meters)", visible=False, placeholder="Enter the length of the area in meters (e.g., 10)")
    area_width = gr.Textbox(label="Width of Area (meters)", visible=False, placeholder="Enter the width of the area in meters (e.g., 5)")
    cell_length = gr.Textbox(label="Cell Opening Length (mm)", visible=False, placeholder="Enter the cell Opening length in millimeters (e.g., 600)")
    cell_width = gr.Textbox(label="Cell Opening Width (mm)", visible=False, placeholder="Enter the cell Opening width in millimeters (e.g., 500)")
    cell_opening_size = gr.Textbox(label="Cell Opening Size (mm)", visible=False, placeholder="Enter the cell opening size in millimeters (e.g., 500)")
    baffle_width = gr.Textbox(label="Baffle Width (mm)", visible=False, placeholder="Enter the baffle width in millimeters (e.g., 25)")

    # Group input fields and buttons together
    with gr.Column():
        grid_type.change(change_inputs, inputs=[grid_type], 
                         outputs=[area_length, area_width, cell_length, cell_width, cell_opening_size, baffle_width])
        with gr.Row():
            submit_btn = gr.Button("Generate Grid Layout", elem_id="generate-btn")
            clear_btn = gr.Button("Clear")

    # Output fields for generated grid (initially hidden)
    output_text = gr.Markdown(visible=False)                              # 1. Text (ok)
    output_image = gr.Image(type="filepath", visible=False)               # 2. Image (ok)
    output_table = gr.HTML(label="Table Preview", visible=False)          # 3. ✅ Must be gr.HTML for html_output
    output_excel = gr.File(label="Download Excel", visible=False)         # 4. Excel path (ok)
    output_file = gr.File(label="Download PDF", visible=False)            # 5. PDF path (ok)
    output_dxf = gr.File(label="Download DXF", visible=False)             # 6. DXF path (ok)
    #table_download = gr.File(label="Download Table (Excel)", visible=False)  # ⬅️ Add this here

    submit_btn.click(
        draw_grid,
        inputs=[grid_type, area_length, area_width, cell_length, cell_width, cell_opening_size, baffle_width],
        outputs=[
            output_text,    # 1
            output_image,   # 2
            output_table,   # 3 ← must be gr.HTML
            output_excel,   # 4
            output_file,    # 5
            output_dxf      # 6
        ]
    )

    # After the button is clicked, we make the output blocks visible
    submit_btn.click(
        lambda: [gr.update(visible=True), gr.update(visible=True), gr.update(visible=True), gr.update(visible=True), gr.update(visible=True), gr.update(visible=True)],
        outputs=[output_text, output_image, output_file, output_dxf, output_table, output_excel]
    )

    # Clear button functionality
    clear_btn.click(
        clear_inputs,
        outputs=[area_length, area_width, cell_length, cell_width, cell_opening_size, baffle_width, output_text, output_image, output_file, output_dxf, output_table,]
    )
    


# Add custom CSS to style the "Generate Grid Layout" button
interface.css = """
#generate-btn {
    background-color: #753bbd !important;
    color: white !important;
    font-weight: bold !important;
    border: none !important;
    box-shadow: none !important;
    outline: none !important;
}
"""

interface.launch()  # Creates a public link