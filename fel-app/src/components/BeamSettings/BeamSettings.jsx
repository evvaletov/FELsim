import ExcelUploadButton from "../ExcelUploadButton/ExcelUploadButton";
import './BeamSettings.css';
import { Container, Row, Button } from "react-bootstrap";

const BeamSettings = ({ setSelectedMenu, excelToAPI, sInterval, setSInterval,
    beamlistSelected, getBeamline }) => {
    return (
        <Container>
            <h4>Graph Settings</h4>
            <ExcelUploadButton excelToAPI={excelToAPI} />
            <label htmlFor="interval" className="forLabels">S axis interval</label>
            <input defaultValue={sInterval}
                    type="number"
                    name="interval" 
                    onChange={(e) => {
                        const val = parseFloat(e.target.value);
                        if (!isNaN(val)) setSInterval(val);
                    }}
            />
            <Row className="mt-2">
                <Button
                    variant="light"
                    onClick={() => {
                        setSelectedMenu(null);
                        getBeamline(beamlistSelected);
                    }}
                >
                    Simulate
                </Button>
            </Row>
        </Container>
    )
};

export default BeamSettings;