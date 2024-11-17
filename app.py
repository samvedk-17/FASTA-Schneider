from flask import Flask , request , redirect , url_for , render_template , send_file
from werkzeug.utils import secure_filename
from Bio import SeqIO
from Bio.Seq import Seq
import os

app = Flask(__name__)

def count_sequences(fasta_file):
    return sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))

    
@app.route("/", methods = ["GET","POST"])
def upload_file():
    if request.method == "POST":
        #Save uploaded file
        file = request.files["file"] 
        if file:
            filepath = secure_filename(file.filename)
            file.save(filepath)

        # Count sequences in the original file
            initial_count = count_sequences(filepath)    

        #Get start and end inputs from frontend
        start = request.form.get("start")
        end = request.form.get("end") 
        length = request.form.get("length")
        length_condition = request.form.get("length_condition")  
        omit_missing = request.form.get("omit_missing") == "on" #checkbox checked status

        #convert length to integer if specified
        length = int(length) if length else None

        #Output file 
        output_path = f"{os.path.splitext(filepath)[0]}_sliced_output.fasta"


        slice_sequences(filepath, output_path, start=start, end=end , length=length, 
                        length_condition=length_condition, omit_missing=omit_missing)

        # Count sequences in the original file
        final_count = count_sequences(output_path)

        preview_lines = []
        with open(output_path , "r") as file:
            for i, line in enumerate(file):
                preview_lines.append(line.strip())
                if i>=13:
                    break
    
        
        #return file for download
        return render_template("results.html" , initial_count=initial_count, final_count=final_count, output_file=output_path , preview=preview_lines)
    
    return render_template("index.html")    
    

        
def slice_sequences(input_fasta , output_fasta , start=None , end = None , length = None, length_condition=None, omit_missing=False):
    sliced_sequences = []

    with open(input_fasta , "r") as input_file:
        for record in SeqIO.parse(input_file , "fasta"):
            sequence = str(record.seq)

            #skip sequences if omit_x box is checked
            if omit_missing and 'X' in sequence:
                continue

            #Determine start and end index  
            start_index = sequence.find(start) if start else 0
            end_index = sequence.find(end , start_index) + len(end) if end else len(sequence)  

            #slice sequence
            if 0 <= start_index < end_index <= len(sequence):
                sliced_seq = sequence[start_index:end_index] 
                #Check length conditions. Code written to skip the seqs not meeting the conditions
                seq_length = len(sliced_seq)
                if (length_condition == "equal" and seq_length != length) or \
                   (length_condition == "greater_than_equal" and seq_length < length) or \
                   (length_condition == "less_than_equal" and seq_length > length):
                    continue
                # Create a new SeqRecord with the sliced sequence
                record.seq = Seq(sliced_seq)
                sliced_sequences.append(record)


    with open(output_fasta , "w") as output_file:
        SeqIO.write(sliced_sequences , output_file , "fasta")

  

@app.route("/download/<path:filepath>")
def download_file(filepath):
    return send_file(filepath, as_attachment=True)

@app.route('/Contact')
def contact():
    return render_template('contact.html')

@app.route('/About')
def about():
    return render_template('about.html')


if __name__ == '__main__':
    app.run(debug=True)

