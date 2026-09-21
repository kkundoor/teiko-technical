PRAGMA foreign_keys = ON;

CREATE TABLE subjects (
    subject TEXT PRIMARY KEY,
    project TEXT NOT NULL,
    condition TEXT NOT NULL,
    age INTEGER NOT NULL CHECK (age >= 0),
    sex TEXT NOT NULL,
    treatment TEXT NOT NULL,
    response TEXT
);

CREATE TABLE samples (
    sample TEXT PRIMARY KEY,
    subject TEXT NOT NULL,
    sample_type TEXT NOT NULL,
    time_from_treatment_start INTEGER NOT NULL,
    FOREIGN KEY (subject)
        REFERENCES subjects(subject)
);

CREATE TABLE cell_counts (
    sample TEXT NOT NULL,
    population TEXT NOT NULL CHECK (
        population IN (
            'b_cell',
            'cd8_t_cell',
            'cd4_t_cell',
            'nk_cell',
            'monocyte'
        )
    ),
    count INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample, population),
    FOREIGN KEY (sample)
        REFERENCES samples(sample)
        ON DELETE CASCADE
);

CREATE INDEX idx_samples_subject
    ON samples(subject);
