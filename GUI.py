from tkinter import *
import tkFileDialog
import tkMessageBox
import os
import fb
import numpy as np
from PIL import ImageTk, Image

MAIN_COLOR = "purple"
FONT = "Callibri"
BASE_COLOR = "white"
UNPRSSED = "light blue"
PRESSED = "pink"
class GUI:
    def __init__(self):
        self.file=""
        self.entry_val=""

        self.main_frame=""
        self.window = Tk()
        self.window.title("Secondary structure prediction")
        self.window.configure(background=BASE_COLOR)
        self.helix_image=self.image_proccess("alpha1.png")
        self.other_image=self.image_proccess("beta.png")
        self.beta_image=self.image_proccess("other (1).PNG")
        self.images_dict={"O":self.other_image,
                     "H":self.helix_image,
                     "E":self.beta_image }


        self.gui_create()
        self.window.mainloop()

    def image_proccess(self, path,width = 15,height = 20 ):
        image = Image.open(path)
        image = image.resize(( height,width), Image.ANTIALIAS) #The (250, 250) is (height, width)
        return ImageTk.PhotoImage(image)

    def run_prdiction(self):
        self.entry_val=self.entry_window.get()
        if(not self.file and not self.entry_val):
            tkMessageBox.showinfo("No input found", "The protein won't predict itself!")
        else:
            seq = self.entry_val if not self.file else self.get_seq_from_file()
            seq=seq.strip("\n\r ")
            if(self.is_seq_legal(seq)):
                self.destroy_frame(self.main_frame)
                self.result_gui(seq)
            else:
                tkMessageBox.showinfo("illegal input", "This is not a legal sequence")

    def is_seq_legal(self,seq):
        '''
        :param seq:
        :return: true if seq is a legal sequnce: aa name, non numeric and upper case
        '''

        not_aa_arr=["B","J","O","U","X","Z"]
        a = [aa for aa in not_aa_arr if aa in seq]


        return seq.isalpha() and seq.isupper() and not a

    def get_seq_from_file(self):
        '''
        called if the user entered a file
        :return: the sequence in that file as a string
        '''
        # f=open(self.file)
        return "".join(self.file.readlines())

    def output_prediction(self,seq,struct):
        images=[]
        row = 2
        col = 1
        for  letter in struct:
            images.append(self.images_dict[letter])
        for image, letter in zip(images,seq):
            aa_label = self.create_label(letter)
            tkImage=Label(self.main_frame, image=image,bd = 0)
            aa_label.grid(row = row , column =col)
            tkImage.grid(row = row+1 , column =col)

            col+=1
            if(col==30):
                col = 1
                aa_label = self.create_label(" ",20)
                aa_label.grid(row = row+2 , columnspan =10)
                row+=3





    def result_gui(self,seq):
        ######Frames######
        self.main_frame = Frame(self.window, bg = "white")
        self.main_frame.grid(row=0)
        ##################
        self.output_prediction(seq,predict(seq))
        empty_line = self.create_label(" ",pady=20)
        empty_line.grid(row =1, columnspan = 1)
        back_B=Button(self.main_frame, text = "return" , command = self.restart ,font=FONT, fg = BASE_COLOR  , activeforeground = BASE_COLOR,bg=UNPRSSED, activebackground=PRESSED)
        back_B.grid(sticky=S,columnspan = 20, pady=70)

    def restart(self):
        self.destroy_frame(self.main_frame)
        self.gui_create()

    def destroy_frame(self,frame):
        for widget in frame.winfo_children():
            widget.destroy()

    def execute(self):
        self.run_prdiction()


    def recive_file(self):

        self.file = tkFileDialog.askopenfile(mode='rb',title='Choose a file')

        if(self.file):

            file_name=self.file.name
            file_name=os.path.basename(file_name)
            received_file = self.create_label("file %s received"%file_name)#Label(self.main_frame,text = "file %s received"%file_name,
                                  #bg = BASE_COLOR,fg= "green",
                                  #font=FONT)
            received_file.grid(row=2,columnspan=2,sticky=W+E+N+S)

    def gui_create(self):
        ######Frames######
        self.main_frame = Frame(self.window, bg = "white")
        self.main_frame.grid(row=0)
        ##################



        #general instruction
        instructions =self.create_label("Browse a file or Enter amino acids sequence:",pady=30,padx=50)# Label(self.main_frame,text = "Browse a file or Enter amino acids sequence:",
        #                      bg = BASE_COLOR ,fg= MAIN_COLOR,
        #                      pady=3,padx = 50 ,height = 3,
        #                      font=FONT)
        or_separator = self.create_label("or",20)#Label(self.main_frame,text = "or",
        #                      bg = BASE_COLOR ,fg= MAIN_COLOR,
        #                      pady=3,height = 3,
        #                      font=FONT)
        new_line1=self.create_label(" ")#Label(self.main_frame,text = " ",bg = BASE_COLOR)
        new_line2=self.create_label(" ")#Label(self.main_frame,text = " ",bg = BASE_COLOR)

        #buttons
        browse_B=Button(self.main_frame, text = "browse file" , command = self.recive_file,font=FONT, fg = BASE_COLOR  , activeforeground = BASE_COLOR,bg=UNPRSSED, activebackground=PRESSED)
        exacute_b = Button(self.main_frame,text = "show me the structure",
                           command = self.execute,
                           bd = 2, height=3, pady=5,relief=RAISED,
                           activeforeground = BASE_COLOR ,activebackground=PRESSED , font=FONT,fg = BASE_COLOR ,bg=UNPRSSED)

        #entry
        input_label =self.create_label("  AA sequence:  ",3)# Label(self.main_frame,text="  AA sequence:  ",bg = BASE_COLOR ,fg= MAIN_COLOR,pady=3,height = 1, font=FONT)
        self.entry_window=Entry(self.main_frame,bg=UNPRSSED )

        #placment
        instructions.grid(row = 0, columnspan = 2)
        browse_B.grid(row = 1, columnspan = 2)
        or_separator.grid(row = 3, columnspan = 2)
        input_label.grid(row = 4, column = 0 )
        self.entry_window.grid(row =4 , column =1)
        new_line1.grid(row = 5, columnspan = 2)
        exacute_b.grid(row = 6, columnspan=2)
        new_line2.grid(row = 7, columnspan = 2)

    def create_label(self,text,pady = 1,padx=1):
        return Label(self.main_frame,text=text,bg = BASE_COLOR ,fg= MAIN_COLOR,pady=pady,height = 1, font=FONT, padx=padx)



def predict(seq):
    e_f = open("emission.txt")
    t_f= open("transition.txt")
    emission = create_nparray(e_f)
    transition = create_nparray(t_f)
    # seq_viterbi = vt.viterbi_round(e, tau, cur_seq)

    return fb.fb_round(emission, transition, seq)[2]

def create_nparray(file):
    matrix=[]
    # file_lines = file.readlines()
    # for line in file_lines:
    #     matrix.append(line.split())
    # return np.array(matrix)
    return np.loadtxt(file)
