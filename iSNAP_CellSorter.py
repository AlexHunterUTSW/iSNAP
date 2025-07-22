import random
import json
import os
from PyQt5.QtWidgets import (
    QWidget, 
    QVBoxLayout, 
    QHBoxLayout,
    QPushButton, 
    QLabel, 
    QGridLayout
)
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QFont

class CellSortGame(QWidget):
    """
    A single-cell sequencing themed minigame implemented using PyQt5.
    Players click on 'diseased' cells as fast as possible.
    The game tracks elapsed time and saves the fastest time locally.
    """
    def __init__(self):
        """
        Initializes the game window and sets up UI components and game state.
        """
        super().__init__()
        self.setWindowTitle("iSNAP Loading... - Cell Sort Challenge")
        self.setGeometry(100, 100, 650, 700) # x, y, width, height
        self.setFixedSize(650, 700) # Prevent resizing for consistent layout

        self.game_active = False
        self.stillLoading = True
        self.time_elapsed = 0.0
        self.game_timer = QTimer(self)
        self.game_timer.timeout.connect(self._update_timer)

        self.cells = [] # Stores references to cell QPushButton objects
        self.target_cells_indices = [] # Stores indices of diseased cells
        self.found_target_cells = 0

        # Game configuration
        self.MAX_GAME_TIME = 60.0 # Maximum game duration in seconds
        self.TOTAL_CELLS = 42 # Total number of cells in the grid
        self.INITIAL_TARGET_CELLS = 5 # Number of diseased cells to find
        self.HIGH_SCORE_FILE = 'cell_sort_high_score.json'
        self.fastest_time = float('inf') # Initialize fastest time to infinity

        self._init_ui()
        self._load_high_score()
        self._display_message("Click 'Start Game' to begin!")

    def _init_ui(self):
        """
        Sets up the main user interface layout and widgets.
        """
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(20, 20, 20, 20)
        main_layout.setSpacing(15)

        # Title
        title_label = QLabel("iSNAP is Loading...")
        title_label.setFont(QFont("Inter", 24, QFont.Bold))
        title_label.setAlignment(Qt.AlignCenter)
        title_label.setStyleSheet("color: #2d3748;") # gray-800
        main_layout.addWidget(title_label)

        # Instructions
        instructions_label = QLabel("Click on all the <span style='color:#ef4444; font-weight:600;'>mutated</span> nucleic acids as fast as you can!")
        instructions_label.setFont(QFont("Inter", 12))
        instructions_label.setAlignment(Qt.AlignCenter)
        instructions_label.setStyleSheet("color: #4a5568;") # gray-600
        main_layout.addWidget(instructions_label)

        # Stats display (Time Elapsed and Fastest Time)
        stats_layout = QHBoxLayout()
        stats_widget = QWidget()
        stats_widget.setLayout(stats_layout)
        stats_widget.setStyleSheet("""
            QWidget {
                background-color: #e2e8f0; /* gray-100 */
                border-radius: 12px; /* rounded-lg */
                padding: 15px;
            }
        """)

        self.timer_label = QLabel("Time Elapsed: <span style='color:#4f46e5;'>0.00</span>s")
        self.timer_label.setFont(QFont("Inter", 14, QFont.DemiBold))
        self.timer_label.setStyleSheet("color: #4a5568;") # gray-700
        stats_layout.addWidget(self.timer_label, alignment=Qt.AlignLeft)

        self.high_score_label = QLabel("Fastest Time: <span style='color:#4f46e5;'>N/A</span>")
        self.high_score_label.setFont(QFont("Inter", 14, QFont.DemiBold))
        self.high_score_label.setStyleSheet("color: #4a5568;") # gray-700
        stats_layout.addWidget(self.high_score_label, alignment=Qt.AlignRight)

        main_layout.addWidget(stats_widget)

        # Game Grid
        self.grid_widget = QWidget()
        self.cell_grid_layout = QGridLayout()
        self.cell_grid_layout.setSpacing(10) # Gap between cells
        self.grid_widget.setLayout(self.cell_grid_layout)
        self.grid_widget.setStyleSheet("""
            QWidget {
                background-color: #e2e8f0; /* Lighter background for the grid area */
                border-radius: 1rem; /* rounded-xl */
                padding: 15px;
                min-height: 250px;
            }
        """)
        main_layout.addWidget(self.grid_widget)

        # Message Box
        self.message_box = QLabel("")
        self.message_box.setFont(QFont("Inter", 11, QFont.Medium))
        self.message_box.setAlignment(Qt.AlignCenter)
        self.message_box.setStyleSheet("""
            QLabel {
                background-color: #dbeafe; /* Light blue */
                color: #1e40af; /* Darker blue */
                padding: 10px;
                border-radius: 8px;
                margin-top: 10px;
            }
        """)
        main_layout.addWidget(self.message_box)

        # Start/Retry Button
        self.start_button = QPushButton("Start Game")
        self.start_button.setFont(QFont("Inter", 13, QFont.DemiBold))
        self.start_button.clicked.connect(self._start_game)
        self.start_button.setStyleSheet("""
            QPushButton {
                background-color: #4f46e5; /* Indigo */
                color: white;
                padding: 12px 24px;
                border-radius: 12px; /* rounded-xl */
                border: none;
            }
            QPushButton:hover {
                background-color: #4338ca; /* Darker indigo */
            }
            QPushButton:pressed {
                background-color: #3730a3; /* Even darker */
            }
            QPushButton:disabled {
                background-color: #94a3b8; /* gray-400 */
                color: #cbd5e1; /* gray-300 */
            }
        """)
        main_layout.addWidget(self.start_button)

        self.setLayout(main_layout)
        self.setStyleSheet("""
            QWidget {
                background-color: #f0f4f8; /* Light blue-gray background */
                border-radius: 1.5rem; /* Rounded corners for the main container */
            }
            CellSortGame {
                background-color: #ffffff; /* White background for the game container itself */
                border-radius: 1.5rem;
            }
        """)

    def _load_high_score(self):
        """
        Loads the fastest time from a local JSON file.
        """
        if os.path.exists(self.HIGH_SCORE_FILE):
            try:
                with open(self.HIGH_SCORE_FILE, 'r') as f:
                    data = json.load(f)
                    self.fastest_time = data.get('fastestTime', float('inf'))
            except json.JSONDecodeError:
                print(f"Error decoding JSON from {self.HIGH_SCORE_FILE}. Resetting high score.")
                self.fastest_time = float('inf')
        self._update_high_score_display()

    def _save_high_score(self):
        """
        Saves the current fastest time to a local JSON file.
        """
        try:
            with open(self.HIGH_SCORE_FILE, 'w') as f:
                json.dump({'fastestTime': self.fastest_time}, f)
        except IOError as e:
            print(f"Error saving high score to {self.HIGH_SCORE_FILE}: {e}")

    def _update_high_score_display(self):
        """
        Updates the high score label in the UI.
        """
        if self.fastest_time == float('inf'):
            self.high_score_label.setText("Fastest Time: <span style='color:#4f46e5;'>N/A</span>")
        else:
            self.high_score_label.setText(f"Fastest Time: <span style='color:#4f46e5;'>{self.fastest_time:.2f}</span>s")

    def _start_game(self):
        """
        Starts a new game round.
        Resets game state, generates new cells, and starts the timer.
        """
        self.game_timer.stop() # Stop any running timer
        self.game_active = True
        self.time_elapsed = 0.0
        self.found_target_cells = 0
        self.cells = []
        self.target_cells_indices = []

        # Clear existing cells from the grid layout
        for i in reversed(range(self.cell_grid_layout.count())):
            widget = self.cell_grid_layout.itemAt(i).widget()
            if widget:
                widget.deleteLater()

        self._display_message("") # Hide any previous messages
        self._generate_cells()
        self.game_timer.start(10) # Start timer updating every 10ms (for 100th of a second precision)

        self.start_button.setText("Restart Game")
        self.start_button.setStyleSheet("""
            QPushButton {
                background-color: #64748b; /* gray-500 */
                color: white;
                padding: 12px 24px;
                border-radius: 12px;
                border: none;
            }
            QPushButton:hover {
                background-color: #475569; /* gray-600 */
            }
            QPushButton:pressed {
                background-color: #334155;
            }
            QPushButton:disabled {
                background-color: #94a3b8;
                color: #cbd5e1;
            }
        """)

    def _generate_cells(self):
        """
        Generates the grid of cells and randomly assigns target (diseased) cells.
        """
        all_cell_indices = list(range(self.TOTAL_CELLS))
        random.shuffle(all_cell_indices)
        self.target_cells_indices = sorted(all_cell_indices[:self.INITIAL_TARGET_CELLS])

        rows = int(self.TOTAL_CELLS**0.5) # Approximate square grid
        cols = (self.TOTAL_CELLS + rows - 1) // rows # Ensure all cells fit

        for i in range(self.TOTAL_CELLS):
            cell_button = QPushButton()
            cell_button.setFixedSize(45, 45) # Cell size
            cell_button.setFont(QFont("Inter", 8))
            cell_button.setProperty("index", i) # Custom property to store index

            is_target = i in self.target_cells_indices
            cell_button.setProperty("is_target", is_target) # Custom property for target status

            cell_button.setStyleSheet(self._get_cell_style(is_target=is_target))
            cell_button.clicked.connect(lambda checked, btn=cell_button: self._handle_cell_click(btn))
            self.cells.append(cell_button)

            row = i // cols
            col = i % cols
            self.cell_grid_layout.addWidget(cell_button, row, col, Qt.AlignCenter)

    def _get_cell_style(self, is_target=False, clicked_correct=False, clicked_incorrect=False):
        """
        Returns the CSS style string for a cell button based on its state.
        """
        base_style = """
            QPushButton {
                border-radius: 22px; /* 50% for circle */
                border: none;
                color: rgba(255, 255, 255, 0.8);
            }
            QPushButton:hover {
            }
        """
        if clicked_correct:
            return base_style + """
                QPushButton { background-color: #22c55e; /* Green */ }
            """
        elif clicked_incorrect:
            return base_style + """
                QPushButton { background-color: #f97316; /* Orange */ }
            """
        elif is_target:
            return base_style + """
                QPushButton { background-color: #ef4444; /* Red */ }
            """
        else:
            return base_style + """
                QPushButton { background-color: #94a3b8; /* Healthy cell color */ }
            """

    def _handle_cell_click(self, clicked_button):
        """
        Handles a click event on a cell button.
        Checks if the clicked cell is a target cell and provides visual feedback.
        """
        if not self.game_active:
            return

        cell_index = clicked_button.property("index")
        is_target = clicked_button.property("is_target")

        # Disable further clicks on this button
        clicked_button.setEnabled(False)

        if is_target:
            self.found_target_cells += 1
            clicked_button.setStyleSheet(self._get_cell_style(clicked_correct=True))
            clicked_button.setText("✓")

            # Check if all target cells have been found
            if self.found_target_cells == self.INITIAL_TARGET_CELLS:
                self._end_game("You found all diseased cells! Well done!", success=True)
        else:
            clicked_button.setStyleSheet(self._get_cell_style(clicked_incorrect=True))
            clicked_button.setText("X")
            # Optional: Add shake animation (requires more advanced PyQt animation)

    def _update_timer(self):
        """
        Updates the game timer every 10 milliseconds.
        Ends the game if the maximum time limit is reached.
        """
        self.time_elapsed += 0.01 # Increment by 0.01 seconds
        self.timer_label.setText(f"Time Elapsed: <span style='color:#4f46e5;'>{self.time_elapsed:.2f}</span>s")

        if self.time_elapsed >= self.MAX_GAME_TIME:
            self._end_game("Time limit reached! Try again.", success=False)

    def _end_game(self, message, success):
        """
        Ends the current game round.
        Stops the timer, displays the final outcome message, and handles high score.
        """
        self.game_timer.stop()
        self.game_active = False

        final_message = f"{message} Your time: {self.time_elapsed:.2f}s"
        self._display_message(final_message)

        if success:
            if self.time_elapsed < self.fastest_time:
                self.fastest_time = self.time_elapsed
                self._save_high_score()
                self._update_high_score_display()
                self._display_message(f"New Fastest Time! {self.fastest_time:.2f}s", is_new_high_score=True)

        # Disable all remaining cells
        for cell_button in self.cells:
            cell_button.setEnabled(False)
            # Ensure target cells that weren't clicked are still visible as red
            if cell_button.property("is_target") and cell_button.text() == "":
                cell_button.setStyleSheet(self._get_cell_style(is_target=True))

        if self.stillLoading:
            self.start_button.setText("Retry Game")
        self.start_button.setStyleSheet("""
            QPushButton {
                background-color: #4f46e5; /* Indigo */
                color: white;
                padding: 12px 24px;
                border-radius: 12px; /* rounded-xl */
                border: none;
            }
            QPushButton:hover {
                background-color: #4338ca; /* Darker indigo */
            }
            QPushButton:pressed {
                background-color: #3730a3; /* Even darker */
            }
        """)

    def _display_message(self, msg, is_new_high_score=False):
        """
        Displays a message in the message box.
        """
        self.message_box.setText(msg)
        if msg:
            self.message_box.show()
            if is_new_high_score:
                self.message_box.setStyleSheet("""
                    QLabel {
                        background-color: #d1fae5; /* Light green */
                        color: #065f46; /* Darker green */
                        padding: 10px;
                        border-radius: 8px;
                        margin-top: 10px;
                    }
                """)
            else:
                self.message_box.setStyleSheet("""
                    QLabel {
                        background-color: #dbeafe; /* Light blue */
                        color: #1e40af; /* Darker blue */
                        padding: 10px;
                        border-radius: 8px;
                        margin-top: 10px;
                    }
                """)
        else:
            self.message_box.hide()
        
    def _load_finished(self):
        if self.game_active:
            self.start_button.setText('Close')
            self.start_button.clicked.connect(self.close)
            self.stillLoading = False
        else:
            self.close()